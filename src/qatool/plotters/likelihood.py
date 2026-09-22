"""Read-support likelihood report (text) over recruited segments.

Kept separate from the plotting code: `report_plot` invokes this only when its
`verbose` flag is set, and callers can run it standalone to get the numbers.
"""

import math
import sqlite3

from qatool.plotters.db import _stripe_coords


def likelihood_report(db_path: str, read_classes: dict,
                      asm1_segs: list, asm2_segs: list,
                      asm1: str, asm2: str,
                      pid_threshold: float = 0.0,
                      plen_threshold: float = 0.0) -> None:
    """
    Compute and print read-support likelihood statistics.

    Pr(r | A) = (1 / L_A) * Σ_p exp(-D(r, p_i))
    where D is the edit distance (NM tag) and L_A is the total length of
    recruited segments for assembly A.

    Likelihood is computed over reads in: asm1_best, asm2_best, common, other.
    asm1_only / asm2_only are counted and reported but excluded from the sum.
    """
    _CLASSES = ["asm1_only", "asm2_only", "asm1_best", "asm2_best", "common", "other"]
    counts = {k: len(read_classes.get(k, set())) for k in _CLASSES}

    L_A1 = sum(s.end - s.start for s in asm1_segs)
    L_A2 = sum(s.end - s.start for s in asm2_segs)

    included = (read_classes.get("asm1_best", set()) |
                read_classes.get("asm2_best", set()) |
                read_classes.get("common",    set()) |
                read_classes.get("other",     set()))

    plen_thr = plen_threshold if plen_threshold is not None else 0.0

    ed_by_asm: dict = {asm1: {}, asm2: {}}
    if included:
        reads_list = list(included)
        con = sqlite3.connect(db_path)
        try:
            for i in range(0, len(reads_list), 800):
                chunk = reads_list[i:i + 800]
                reads_ph = ",".join(["?"] * len(chunk))
                sql = f"""
                    SELECT assembly, read_name, ed FROM aln
                    WHERE assembly IN (?, ?)
                      AND contig IS NOT NULL
                      AND pid >= ? AND plen >= ?
                      AND ed IS NOT NULL
                      AND read_name IN ({reads_ph})
                """
                for row_asm, rname, ed in con.execute(
                        sql, [asm1, asm2, pid_threshold, plen_thr] + chunk):
                    if row_asm in ed_by_asm:
                        ed_by_asm[row_asm].setdefault(rname, []).append(int(ed))
        finally:
            con.close()

    def _logsumexp_neg(eds):
        m = min(eds)
        return -m + math.log(sum(math.exp(-(e - m)) for e in eds))

    def _total_log_likelihood(asm_ed_map, L_A):
        if L_A <= 0:
            return float('-inf'), len(included)
        ll, n_missing = 0.0, 0
        log_L_A = math.log(L_A)
        for rname in included:
            eds = asm_ed_map.get(rname)
            if not eds:
                n_missing += 1
                continue
            ll += _logsumexp_neg(eds) - log_L_A
        return ll, n_missing

    ll1, miss1 = _total_log_likelihood(ed_by_asm[asm1], L_A1)
    ll2, miss2 = _total_log_likelihood(ed_by_asm[asm2], L_A2)
    n_inc = len(included)
    delta = ll2 - ll1 if math.isfinite(ll1) and math.isfinite(ll2) else float('nan')

    w = max(len(asm1), len(asm2), 12)
    print()
    print("=== Likelihood Report ===")
    print("Read class counts:")
    for cls in _CLASSES:
        print(f"  {cls:<12} : {counts[cls]:>7,}")
    print(f"  {'total':<12} : {sum(counts.values()):>7,}")
    print()
    print("Recruited segment lengths:")
    print(f"  {asm1:<{w}} : {L_A1:>14,} bp")
    print(f"  {asm2:<{w}} : {L_A2:>14,} bp")
    print()
    def _fmt_coords(coord_map):
        parts = [f"{c}:{s}-{e}({n})" for c, ivs in coord_map.items() for s, e, n in ivs]
        return "[" + ", ".join(parts) + "]" if parts else "[]"

    asm1_coords = _stripe_coords(db_path, read_classes.get("asm1_only", set()), asm1, pid_threshold, plen_threshold)
    asm2_coords = _stripe_coords(db_path, read_classes.get("asm2_only", set()), asm2, pid_threshold, plen_threshold)
    print(f"{asm1} support only: {_fmt_coords(asm1_coords)}")
    print(f"{asm2} support only: {_fmt_coords(asm2_coords)}")
    print()
    print(f"Log-likelihood  (n={n_inc:,} reads: asm1_best + asm2_best + common + other)")
    print(f"  log P(reads | {asm1:<{w}}) = {ll1:>14.3f}   [reads without aln: {miss1}]")
    print(f"  log P(reads | {asm2:<{w}}) = {ll2:>14.3f}   [reads without aln: {miss2}]")
    if math.isfinite(delta):
        better = asm2 if delta > 0 else asm1
        sign = "+" if delta > 0 else ""
        print(f"  ΔLLR  ({asm2} − {asm1}){' ' * (w - len(asm2))} = {sign}{delta:>13.3f}   [{better} better supported]")
    else:
        print(f"  ΔLLR  ({asm2} − {asm1}) = N/A")
    print()


def phmm_likelihood_report(segments, read_names, asm1: str, asm2: str, *,
                           db_path: str, bams: dict, fastas: dict,
                           genome_fasta: str, library_assembly: str,
                           C: float = 200.0, **kwargs):
    """Pair-HMM locus score: emission, exposure and Delta S, printed and returned.

    The replacement for `likelihood_report` above, which scores reads by
    exp(-NM) and keeps -n log L_A -- the conditional form that Proposition 3 of
    model.tex shows is blind to copy number. This one sums the pair HMM over
    every window a read has in each candidate, clamps each read at +/-C, and adds
    the exposure term -2 lambda (L_eff_1 - L_eff_2).

    `genome_fasta` and `library_assembly` fix lambda: read starts per base per
    strand, from the primary reads aligned to `library_assembly` over the diploid
    length of `genome_fasta`. Use the same pair for every locus you compare --
    lambda must not depend on the locus under test.
    """
    from qamodel.exposure import fasta_length, library_stats
    from qamodel.score import score_locus

    lib = library_stats(db_path, library_assembly, fasta_length(genome_fasta))
    result = score_locus(read_names, segments, asm1, asm2, db_path=db_path,
                         bams=bams, fastas=fastas, library=lib, C=C, **kwargs)
    print(result.summary())
    return result
