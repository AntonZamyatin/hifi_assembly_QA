"""Tests for the reads under the alignment strip (qaview.pileup).

The synthetic cases are the one that matters on screen: two assemblies that
differ by an AT inside an (AT)n repeat, with reads whose aligner put the indel
somewhere else in the repeat. Every read must show it in the alignment's gap
columns, base for base. The last tests check the projection against real
records and skip when the data is absent.

Run: venv/bin/python tests/test_pileup.py
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam

from qaview.pileup import column_maps, events, project, revcomp

_failures = []

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "examples")


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _record(cigar, start, reverse=False):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "c", "LN": 100000}]})
    a = pysam.AlignedSegment(header)
    a.query_name = "r"
    a.reference_id = 0
    a.reference_start = start
    a.cigarstring = cigar
    a.flag = 16 if reverse else 0
    return a


def states(runs):
    """Runs -> one letter per column, from the first run's column."""
    return "".join(kind * (c1 - c0 + 1) for c0, c1, kind in runs)


# top lacks one AT that bottom has, inside the repeat: the alignment's gap is left-aligned
A1 = "GGCACAT" + "AT" * 4 + "CCGA"          # 19 bp
A2 = "GGCACAT" + "AT" * 5 + "CCGA"          # 21 bp
TOP = "GGCACAT--" + "AT" * 4 + "CCGA"
BOT = A2


def test_column_maps_skip_gaps_and_outside():
    a = {"top": "  AC-GT", "bottom": "TTACGGT", "strand": "+"}
    top, bottom = column_maps(a)
    check(top == [2, 3, 5, 6], f"top bases sit in columns {top}")
    check(bottom == [0, 1, 2, 3, 4, 5, 6], f"bottom bases {bottom}")


def test_insertion_fills_the_gap_columns():
    cols, _ = column_maps({"top": TOP, "bottom": BOT})
    # a read on asm1 carrying the extra AT, aligner put it at the repeat's end
    read = A1[:15] + "AT" + A1[15:]
    ev = events(_record("15M2I4M", 0), read, A1, 0, len(A1))
    runs, ticks = project(ev, cols, BOT)
    check(states(runs) == "b" * 7 + "ii" + "b" * 12 and not ticks,
          f"the inserted AT fills gap columns 7-8: {states(runs)} {ticks}")
    # a read on asm1 exactly as asm1: the gap columns are gaps in the read too
    runs, _ = project(events(_record("19M", 0), A1, A1, 0, len(A1)), cols, BOT)
    check(states(runs) == "b" * 7 + "gg" + "b" * 12, f"no insertion: {states(runs)}")


def test_deletion_lands_on_the_gap_columns():
    _, cols = column_maps({"top": TOP, "bottom": BOT})
    # a read on asm2 lacking one AT, aligner put the deletion at the repeat's end
    read = A2[:15] + A2[17:]
    runs, _ = project(events(_record("15M2D4M", 0), read, A2, 0, len(A2)), cols, TOP)
    check(states(runs) == "b" * 7 + "dd" + "b" * 12, f"deletion under the gap: {states(runs)}")


def test_reverse_complemented_row():
    """The bottom row shown reverse-complemented; reads aligned in genome orientation."""
    g2 = revcomp(A2)                        # the genome; the row shows revcomp(g2) == A2
    _, cols = column_maps({"top": TOP, "bottom": A2})
    # in genome orientation the missing AT is one of the repeat's; put it at the left end
    rep = g2.index("ATAT")
    read = g2[:rep] + g2[rep + 2:]
    runs, _ = project(events(_record(f"{rep}M2D{len(g2) - rep - 2}M", 0), read, g2, 0, len(g2),
                             reverse=True), cols, TOP)
    check(states(runs) == "b" * 7 + "dd" + "b" * 12, f"reverse-complemented: {states(runs)}")
    # and a reverse-strand read carrying one AT fewer, hard-clipped in front
    rec = _record(f"3H{rep}M2D{len(g2) - rep - 2}M", 0, reverse=True)
    as_sequenced = revcomp("NNN" + read)[:len(read)]    # clipped bases are not in the alignment
    runs, _ = project(events(rec, revcomp("TTT" + read), g2, 0, len(g2), reverse=True), cols, TOP)
    check(states(runs) == "b" * 7 + "dd" + "b" * 12, f"reverse strand, hard clip: {states(runs)}")


def test_mismatch_and_novel_insertion():
    cols, _ = column_maps({"top": TOP, "bottom": BOT})
    read = A1[:2] + "T" + A1[3:17] + "GG" + A1[17:]    # a mismatch at 2, GG inserted near the end
    runs, ticks = project(events(_record("17M2I2M", 0), read, A1, 0, len(A1)), cols, BOT)
    st = states(runs)
    check(st[2] == "x" and st.count("x") == 1, f"one mismatch at column 2: {st}")
    check(ticks == [[cols[17], 2]], f"an insertion neither assembly has is a tick: {ticks}")


def _reads(bam, contig, lo, hi, store, limit):
    n = 0
    for a in bam.fetch(contig, lo, hi):
        if n >= limit or not a.has_tag("NM"):
            continue
        try:
            seq, _ = store.get(a.query_name)
        except (KeyError, ValueError):
            continue
        n += 1
        yield a, seq


def _data():
    p = {"bam": os.path.join(DATA, "alignments/hg002.curated.sorted.bam"),
         "fa": os.path.join(DATA, "bgz/hg002.curated.fasta.bgz"),
         "db": os.path.join(DATA, "dbs/hg002.aln.sqlite")}
    return p if all(os.path.exists(v) for v in p.values()) else None


def test_projection_agrees_with_nm_on_real_records():
    """Projected on its own reference with no gap columns, a record's differences are NM."""
    p = _data()
    if p is None:
        print("  SKIP: data absent")
        return
    from qamodel.phmm.reads import ReadStore

    bam, fa = pysam.AlignmentFile(p["bam"]), pysam.FastaFile(p["fa"])
    n = bad = 0
    with ReadStore(p["db"], {"curated": p["bam"]}) as store:
        # chrM: many supplementary, hard-clipped records on both strands
        for a, seq in _reads(bam, "CP139510.1", 0, 16570, store, 150):
            s, e = a.reference_start, a.reference_end
            ref = fa.fetch("CP139510.1", s, e).upper()
            for rev in (False, True):
                runs, ticks = project(events(a, seq, ref, s, e, rev), list(range(e - s)), ref)
                st = states(runs)
                edits = st.count("x") + st.count("d") + sum(t[1] for t in ticks)
                n += 1
                bad += edits != a.get_tag("NM")
    check(n > 200 and bad == 0, f"{bad} of {n} projections disagree with NM")


def test_every_read_base_is_placed_once():
    """Inside a stretch with gap columns, bases shown + ticks == the read's bases there."""
    p = _data()
    if p is None:
        print("  SKIP: data absent")
        return
    from qamodel.phmm.reads import ReadStore

    bam, fa = pysam.AlignmentFile(p["bam"]), pysam.FastaFile(p["fa"])
    lo, hi, contig = 93_000, 95_000, "CP139528.2"
    ref = fa.fetch(contig, lo, hi).upper()
    rng = random.Random(1)
    row = "".join(c + ("-" * 2 if rng.random() < 0.01 else "") for c in ref)   # sprinkle gap columns
    cols = [k for k, c in enumerate(row) if c != "-"]
    other = "".join("A" if c == "-" else c for c in row)
    n = bad = 0
    with ReadStore(p["db"], {"curated": p["bam"]}) as store:
        for a, seq in _reads(bam, contig, lo, hi, store, 80):
            s, e = max(a.reference_start, lo), min(a.reference_end, hi)
            expect = 0
            r = a.reference_start
            for op, ln in a.cigartuples:
                if op in (0, 7, 8):
                    expect += max(0, min(r + ln, e) - max(r, s))
                    r += ln
                elif op == 1:
                    expect += ln if s < r < e else 0
                elif op in (2, 3):
                    r += ln
            runs, ticks = project(events(a, seq, ref, lo, hi), cols, other)
            st = states(runs)
            shown = st.count("b") + st.count("x") + st.count("i") + sum(t[1] for t in ticks)
            n += 1
            bad += shown != expect
    check(n > 20 and bad == 0, f"{bad} of {n} reads misplace bases")


def main():
    tests = [test_column_maps_skip_gaps_and_outside, test_insertion_fills_the_gap_columns,
             test_deletion_lands_on_the_gap_columns, test_reverse_complemented_row,
             test_mismatch_and_novel_insertion, test_projection_agrees_with_nm_on_real_records,
             test_every_read_base_is_placed_once]
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
