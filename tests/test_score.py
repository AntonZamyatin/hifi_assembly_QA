"""Tests for the locus score: clamp, exposure, and the whole thing end to end.

The end-to-end test is Proposition 6 on synthetic sequence. The truth has one
copy of a region, the other candidate falsely duplicates it, and reads are drawn
from the truth. The emission term must come out at -log 2 per read -- favouring
the fabrication -- and the exposure term must be what turns Delta S around. That
is the case a number built from emissions alone gets confidently wrong.

Run: venv/bin/python tests/test_score.py
"""

import math
import os
import sys
from types import SimpleNamespace as NS

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.exposure import Library, effective_length, exposure_term
from qamodel.phmm import PhmmParams, quality_table, sample_read
from qamodel.phmm.kernel import forward_numba
from qamodel.phmm.windows import Placement, build_windows, straddles
from qamodel.score import LocusScore, ReadScore, placements, score_read

QT = quality_table()
_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def point_library(n_reads, genome, length, bin_bp=100):
    """Every read the same length: makes L_eff exact and checkable by hand."""
    hist = np.zeros(length // bin_bp + 1, dtype=np.int64)
    hist[length // bin_bp] = n_reads
    return Library(n_reads, genome, bin_bp, hist)


def test_clamp_and_infinities():
    r = ReadScore("r", 100, -5.0, -math.inf, 1, 0, 0.0)
    check(r.ratio == math.inf and r.clamped(200) == 200, "one-sided read must land on +C")
    r = ReadScore("r", 100, -math.inf, -5.0, 0, 1, 0.0)
    check(r.clamped(200) == -200, "one-sided the other way must land on -C")
    r = ReadScore("r", 100, -5.0, -9.0, 1, 1, 0.0)
    check(r.clamped(200) == 4.0 and r.clamped(3) == 3.0, "clamp must be continuous")


def test_sweep_needs_no_rescoring():
    lib = point_library(10, 10_000, 1000)
    x = exposure_term([NS(start=0, end=5000)], [NS(start=0, end=5000)], lib)
    reads = [ReadScore(str(i), 100, 0.0, -v, 1, 1, 0.0) for i, v in enumerate((1, 50, 500))]
    ls = LocusScore("a", "b", reads, x, C=100.0)
    check(abs(ls.emission - (1 + 50 + 100)) < 1e-12, f"C=100 emission {ls.emission}")
    check(abs(ls.at(1000) - (1 + 50 + 500)) < 1e-12, f"C=1000 gives {ls.at(1000)}")
    check(ls.clamped_counts() == (1, 0), f"{ls.clamped_counts()}")


def test_effective_length_exact_form():
    lib = point_library(100, 10**6, 15_000)
    got = effective_length([45_000], lib)
    want = 45_000 - lib.centres[lib.hist.argmax()] + 1
    check(abs(got - want) < 1e-6, f"L_eff {got} vs L - l + 1 = {want}")
    check(effective_length([10_000], lib) == 0.0,
          "a segment shorter than every read admits no starts")
    two = effective_length([45_000, 45_000], lib)
    check(abs(two - 2 * want) < 1e-6, "fragmentation: two segments pay the read length twice")


def test_lambda_definition():
    lib = point_library(12_561_909, 5_999_424_718, 15_400)
    check(abs(lib.lam - 12_561_909 / (2 * 5_999_424_718)) < 1e-15, "lambda = N / 2 L_dip")
    check(abs(lib.coverage - 2 * lib.lam * lib.mean_length) < 1e-12, "coverage identity")


def test_false_duplication_end_to_end():
    """Proposition 6: emission prefers the fabrication; exposure must reverse it."""
    rng = np.random.default_rng(6)
    R = lambda n: "".join(rng.choice(list("ACGT"), size=n))
    region, fl, fr = R(3000), R(400), R(400)
    A1 = fl + region + fr                        # truth: one copy
    A2 = fl + region + region + fr               # candidate: falsely duplicated
    fa = {"A1": NS(fetch=lambda c, lo, hi: A1[lo:hi]),
          "A2": NS(fetch=lambda c, lo, hi: A2[lo:hi])}
    seg1 = [NS(contig="c", start=0, end=len(A1))]
    seg2 = [NS(contig="c", start=0, end=len(A2))]
    err = PhmmParams(delta_i=2e-3, delta_d=2e-3)
    L = 800

    # Origins stop 100 bp short of the region's end. A read that runs past the
    # duplicated region into the flank is a junction read (Remark 12): only one
    # copy in the fabricated candidate has that flank beside it, so it sees one
    # good placement in each candidate and scores 0, not -log 2. That is the model
    # being right, but it is not the case this test isolates.
    scores = []
    for _ in range(60):
        o = int(rng.integers(400, 400 + 3000 - L - 100))
        s = sample_read(A1, err, rng, origin=o, qual=40, max_len=L)
        if s.truncated:
            continue
        p1 = [Placement("c", o, o + L, 1, f"{L}M", 0, 0)]
        p2 = [Placement("c", o, o + L, 1, f"{L}M", 0, 0),
              Placement("c", o + 3000, o + 3000 + L, 1, f"{L}M", 0, 0)]
        w1 = build_windows(p1, seg1, len(s.read), fa["A1"]).windows
        w2 = build_windows(p2, seg2, len(s.read), fa["A2"]).windows
        lp1, n1, e1 = score_read(s.read, s.qual, w1, err, QT, forward_numba)
        lp2, n2, e2 = score_read(s.read, s.qual, w2, err, QT, forward_numba)
        scores.append(ReadScore(str(_), L, lp1, lp2, n1, n2, max(e1, e2)))

    n = len(scores)
    # the reads the true region draws: n = 2 * lambda * L_eff(A1)
    lib = point_library(1, 1, L)
    leff1 = effective_length([len(A1)], lib)
    lam = n / (2 * leff1)
    lib = Library(n, int(round(n / (2 * lam))), lib.bin_bp, lib.hist)
    ls = LocusScore("truth", "dup", scores, exposure_term(seg1, seg2, lib), C=200.0)

    per_read = ls.emission / n
    check(abs(per_read + math.log(2)) < 0.01,
          f"emission {per_read:+.4f}/read, expected -log 2 = {-math.log(2):+.4f}")
    check(ls.exposure.value > 0, "exposure must favour the shorter, true candidate")
    check(ls.delta_s > 0,
          f"Delta S {ls.delta_s:+.2f} must favour the truth once exposure is in")
    print(f"    {n} reads: emission {ls.emission:+.2f} ({per_read:+.4f}/read), "
          f"exposure {ls.exposure.value:+.2f}, Delta S {ls.delta_s:+.2f}")


def test_boundary_read_scores_as_if_the_unit_held_the_flank():
    """A boundary read crosses a segment edge into flank both candidates share.
    Scored with windows past the edge, its ratio is the one it would have if the
    segments simply contained it -- the flank cancels."""
    rng = np.random.default_rng(12)
    R = lambda n: "".join(rng.choice(list("ACGT"), size=n))
    left, inside, flank = R(500), R(1500), R(1500)
    A1 = left + inside + flank
    A2 = left + inside[:1200] + "TTT" + inside[1200:] + flank  # differs inside the segment only
    fa = {"A1": NS(fetch=lambda c, lo, hi: A1[lo:hi]), "A2": NS(fetch=lambda c, lo, hi: A2[lo:hi])}
    edge1, edge2 = len(left) + len(inside), len(left) + len(inside) + 3
    seg = {"A1": [NS(contig="c", start=0, end=edge1)], "A2": [NS(contig="c", start=0, end=edge2)]}
    whole = {"A1": [NS(contig="c", start=0, end=len(A1))], "A2": [NS(contig="c", start=0, end=len(A2))]}
    err = PhmmParams(delta_i=2e-3, delta_d=2e-3)
    L, o = 1200, 1500                       # 1,500-2,700 in A1: over the edge at 2,000
    sr = sample_read(A1, err, rng, origin=o, qual=40, max_len=L)
    pl = {"A1": [Placement("c", o, o + L, 1, f"{L}M", 0, 0)],
          "A2": [Placement("c", o, o + L + 3, 1, f"{L}M", 0, 0)]}

    def ratio(segs, cross):
        lp = {}
        for a in ("A1", "A2"):
            ws = build_windows(pl[a], segs[a], len(sr.read), fa[a], cross_edges=cross)
            lp[a] = score_read(sr.read, sr.qual, ws.windows, err, QT, forward_numba)[0]
        return lp["A1"] - lp["A2"]

    check(straddles(pl["A1"], seg["A1"]), "the read crosses the edge: an ordinary read is dropped")
    contained, boundary = ratio(whole, False), ratio(seg, True)
    check(abs(contained - boundary) < 1e-9,
          f"boundary ratio {boundary:.6f} equals the contained one {contained:.6f}")
    check(boundary > 5, f"and it sees the 3 bp difference inside: {boundary:.2f} nats for A1")


def test_placements_follow_recruitment_thresholds():
    """Only alignments the recruiter would follow are placements: a junk secondary
    crossing an edge can no longer drop a read whose real placement is inside."""
    import sqlite3
    import tempfile

    from fixtures import write_db

    with tempfile.TemporaryDirectory() as tmp:
        db = sqlite3.connect(write_db(os.path.join(tmp, "a.sqlite"), [
            {"assembly": "A", "read_name": "r", "contig": "c", "start": 1000, "end": 3000,
             "pid": 0.9995, "plen": 1.0},
            {"assembly": "A", "read_name": "r", "contig": "c", "start": 9000, "end": 10_800,
             "pid": 0.95, "plen": 0.9, "aln_type": "secondary"},
        ]))
        seg = [NS(contig="c", start=0, end=10_000)]
        every, _ = placements(db, "A", "r")
        check(len(every) == 2 and straddles(every, seg),
              "without thresholds the junk secondary crosses the edge (the old drop)")
        every, kept = placements(db, "A", "r", pid_min=0.997, plen_min=0.999)
        check([(p.start, p.end) for p in every] == [(1000, 3000)] and not straddles(every, seg),
              "with the recruitment's thresholds only the real placement counts")
        none, _ = placements(db, "A", "r", pid_min=0.9999)
        check(none == [], "nothing passing: no placement, so no window in that candidate")


def test_reads_without_qualities_score_at_a_fixed_quality():
    """Reads aligned from FASTA carry no QUAL: each base gets FALLBACK_QUAL, the read
    is returned as sequenced on either strand, scored, and counted as such."""
    import random
    import tempfile

    import pysam

    from fixtures import write_db
    from qamodel.phmm.reads import FALLBACK_QUAL, ReadStore, revcomp
    from qamodel.score import score_locus

    rng = random.Random(2)
    ref = "".join(rng.choice("ACGT") for _ in range(6_000))
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "g.fa")
        with open(fa, "w") as fh:
            fh.write(f">c\n{ref}\n")
        pysam.faidx(fa)
        bam = os.path.join(tmp, "g.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "c", "LN": len(ref)}]}
        rows = []
        with pysam.AlignmentFile(bam, "wb", header=header) as out:
            for name, start, rev in (("fwd", 1_000, False), ("rev", 3_000, True)):
                a = pysam.AlignedSegment()
                a.query_name, a.reference_id, a.reference_start = name, 0, start
                a.mapping_quality, a.cigartuples, a.flag = 60, [(0, 1_500)], 16 if rev else 0
                a.query_sequence = ref[start:start + 1_500]       # stored on the reference strand
                a.set_tag("NM", 0)
                out.write(a)                                      # no qualities: QUAL '*'
                for asm in ("A", "B"):
                    rows.append({"assembly": asm, "read_name": name, "contig": "c", "start": start,
                                 "end": start + 1_500, "strand": 0 if rev else 1, "q_len": 1_500,
                                 "q_aln_len": 1_500})
        pysam.index(bam)
        db = write_db(os.path.join(tmp, "a.sqlite"), rows)

        with ReadStore(db, {"A": bam}) as store:
            seq, qual = store.get("rev")
            check(seq == revcomp(ref[3_000:4_500]) and (qual == FALLBACK_QUAL).all() and len(qual) == 1_500,
                  "a reverse read comes back as sequenced, every base at the fixed quality")
            store.get("fwd")
            check(store.fixed_quality == {"fwd", "rev"}, f"both noted: {store.fixed_quality}")
        with ReadStore(db, {"A": bam}, fallback_qual=None) as store:
            try:
                store.get("fwd")
                check(False, "without a fallback a read with no qualities is refused")
            except ValueError as exc:
                check("no base qualities" in str(exc), str(exc))

        segs = [NS(assembly=a, contig="c", start=0, end=len(ref)) for a in ("A", "B")]
        res = score_locus(["fwd", "rev"], segs, "A", "B", db_path=db, bams={"A": bam, "B": bam},
                          fastas={"A": fa, "B": fa}, library=point_library(10_000, 1_000_000, 1_500),
                          workers=1)
        check(len(res.reads) == 2 and sorted(res.fixed_quality) == ["fwd", "rev"],
              f"scored, and counted on the score: {res.fixed_quality}")
        check(all(abs(r.ratio) < 1e-9 for r in res.reads), "identical candidates: ratio 0")


def main():
    tests = [test_clamp_and_infinities, test_sweep_needs_no_rescoring,
             test_effective_length_exact_form, test_lambda_definition,
             test_false_duplication_end_to_end,
             test_boundary_read_scores_as_if_the_unit_held_the_flank,
             test_placements_follow_recruitment_thresholds,
             test_reads_without_qualities_score_at_a_fixed_quality]
    for fn in tests:
        print(f"{fn.__name__} ...")
        fn()
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
