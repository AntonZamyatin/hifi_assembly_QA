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


# Bump when a change to the scoring code alters results, so cached scores from
# before it are not reused.
SCORE_CACHE_VERSION = 1


def locus_id(segments, read_names, asm1: str, asm2: str) -> str:
    """Content id of a recruited locus: its segments and read set, not the seed.

    Seeds inside one repeat family usually recruit the same locus; they share
    this id, so a batch can report the locus once and score it once.
    """
    import hashlib
    import json

    segs = sorted((s.assembly, s.contig, int(s.start), int(s.end)) for s in segments)
    payload = json.dumps({"asm": [asm1, asm2], "segments": segs,
                          "reads": sorted(read_names)})
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def _score_key(lid: str, asm1: str, asm2: str, db_path: str, bams: dict, fastas: dict,
               genome_fasta: str, library_assembly: str, C: float, kwargs: dict) -> str:
    """locus_id plus everything else that changes a LocusScore."""
    import hashlib
    import inspect
    import json
    from dataclasses import asdict

    from qamodel.phmm import PhmmParams
    from qamodel.score import score_locus

    sig = inspect.signature(score_locus).parameters
    settings = {n: kwargs.get(n, sig[n].default) for n in ("k", "extra", "width", "prune")}
    payload = json.dumps({
        "version": SCORE_CACHE_VERSION, "locus": lid, "asm": [asm1, asm2],
        "db": db_path, "bams": [bams[asm1], bams[asm2]], "fastas": [fastas[asm1], fastas[asm2]],
        "genome_fasta": genome_fasta, "library_assembly": library_assembly, "C": C,
        "phmm": asdict(kwargs.get("params") or PhmmParams()),
        "qt": repr(kwargs.get("qt")), **settings,
    }, sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def score_cache_path(cache_dir: str, lid: str, asm1: str, asm2: str, *, db_path: str,
                     bams: dict, fastas: dict, genome_fasta: str, library_assembly: str,
                     C: float = 200.0, **kwargs) -> str:
    """Where phmm_likelihood_report(cache_dir=...) keeps the score of locus `lid`.

    Lets a caller see whether a locus is already scored without scoring it.
    """
    import os

    key = _score_key(lid, asm1, asm2, db_path, bams, fastas, genome_fasta,
                     library_assembly, C, kwargs)
    return os.path.join(cache_dir, f"score_{key}.pkl")


def phmm_likelihood_report(segments, read_names, asm1: str, asm2: str, *,
                           db_path: str, bams: dict, fastas: dict,
                           genome_fasta: str, library_assembly: str,
                           C: float = 200.0, cache_dir: str | None = None, **kwargs):
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

    With `cache_dir`, the result is stored under a key made from the locus
    content (`locus_id`) and every setting that changes the score, and reused
    when another seed recruits the same locus. A lock serialises concurrent runs
    on one locus, so a second notebook waits for the first instead of repeating it.
    """
    from qamodel.exposure import fasta_length, library_stats
    from qamodel.score import score_locus

    def compute():
        lib = library_stats(db_path, library_assembly, fasta_length(genome_fasta))
        return score_locus(read_names, segments, asm1, asm2, db_path=db_path,
                           bams=bams, fastas=fastas, library=lib, C=C, **kwargs)

    if cache_dir is None:
        result = compute()
    else:
        import fcntl
        import os
        import pickle

        lid = locus_id(segments, read_names, asm1, asm2)
        path = score_cache_path(cache_dir, lid, asm1, asm2, db_path=db_path, bams=bams,
                                fastas=fastas, genome_fasta=genome_fasta,
                                library_assembly=library_assembly, C=C, **kwargs)
        os.makedirs(cache_dir, exist_ok=True)
        with open(path + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if os.path.exists(path):
                with open(path, "rb") as fh:
                    result = pickle.load(fh)
                print(f"[cache] locus {lid} already scored: {os.path.basename(path)}")
            else:
                result = compute()
                tmp = f"{path}.{os.getpid()}.tmp"
                with open(tmp, "wb") as fh:
                    pickle.dump(result, fh)
                os.replace(tmp, path)
    print(result.summary())
    return result
