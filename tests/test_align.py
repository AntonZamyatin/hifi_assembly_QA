"""Tests for the shift+drag stretch alignment (qaview.align).

Synthetic sequences only: a stub stands in for RegionData (segment layout and
binsize) and for the FASTA files.

Run: venv/bin/python tests/test_align.py
"""

import os
import random
import sys
import time
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import edlib

from qaview.align import MAX_SPAN, _expand, align_rect, layout, left_align, polish, revcomp, span_at

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


class Fasta:
    def __init__(self, seqs):
        self.seqs = seqs

    def fetch(self, contig, lo, hi):
        return self.seqs[contig][lo:hi]


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def plant(seq, deletion_at, deletion_len, substitutions):
    """seq with a deletion and substitutions (positions in the edited sequence)."""
    out = seq[:deletion_at] + seq[deletion_at + deletion_len:]
    for p in substitutions:
        out = out[:p] + ("A" if out[p] != "A" else "C") + out[p + 1:]
    return out


def _setup(asm1, asm2, binsize):
    rd = SimpleNamespace(binsize=binsize, asm1="A", asm2="B",
                         asm1_segs_info=[("a", 0, len(asm1), 0)],
                         asm2_segs_info=[("b", 0, len(asm2), 0)])
    return rd, Fasta({"a": asm1}), Fasta({"b": asm2})


def test_layout_rows_and_counts():
    q, t = "ACGTAACGTTTC", "ACGTACGTTAC"          # one base only in q, one mismatch
    ops = _expand(edlib.align(q, t, mode="NW", task="path")["cigar"])
    a = layout(ops, q, t, 0, 0, len(t))
    check(len(a["top"]) == len(a["mid"]) == len(a["bottom"]), "rows have equal length")
    check(a["top"].replace("-", "") == q and a["bottom"].replace("-", "") == t,
          "rows spell the two sequences")
    check(a["ins"] == 1 and a["del"] == 0 and a["mismatches"] == 1,
          f"counts: {a['ins']} ins, {a['del']} del, {a['mismatches']} mm")


def test_polish_joins_a_scattered_gap():
    """Unit costs may split one 12 bp indel into pieces; polishing makes it one."""
    rng = random.Random(3)
    q = rand_seq(400, rng)
    t = q[:200] + q[212:]
    ops = polish(_expand(edlib.align(q, t, mode="NW", task="path")["cigar"]), q, t)
    runs = "".join(ops).replace("=", " ").split()
    check(runs == ["I" * 12], f"one 12 bp insertion, got runs {runs}")


def test_gaps_are_left_aligned():
    t, m, b = left_align("CATATAT--ATAG", "|||||||  ||||", "CATATATATATAG")
    check((t, m, b) == ("C--ATATATATAG", "|  ||||||||||", "CATATATATATAG"),
          f"a gap in a repeat moves to its left end: {t!r}")
    t, m, b = left_align("CAT--G", "|||  |", "CATCGG")
    check(t == "CAT--G", f"a gap stops where the sequence differs: {t!r}")


def test_span_at_snaps_to_bins_and_keeps_one_segment():
    segs = [("c1", 1000, 1250, 0), ("c2", 500, 700, 250)]   # (contig, start, end, offset)
    check(span_at(segs, 100, 30, 150) == ("c1", 1000, 1250, 1000, 1200, 1090, False),
          "snapped outward to whole bins")
    check(span_at(segs, 100, 170, 170)[3:5] == (1100, 1200), "a point is its bin")
    check(span_at(segs, 100, 240, 245)[3:5] == (1200, 1250), "the last bin is clipped to the segment")
    got = span_at(segs, 100, 100, 330)                      # midpoint 215 lies in c1
    check(got[0] == "c1" and got[3:] == (1100, 1250, 1215, True),
          f"a rectangle across segments keeps the one under its midpoint: {got}")


def test_whole_stretches_with_planted_variants():
    """Both whole stretches are shown along their true offset; exactly the variants count."""
    rng = random.Random(4)
    asm1 = rand_seq(8000, rng)
    asm2 = rand_seq(330, rng) + plant(asm1, 3500, 12, (2600, 4100))
    rd, f1, f2 = _setup(asm1, asm2, binsize=1000)
    a = align_rect(rd, f1, f2, x0=2400, x1=4900, y0=2100, y1=4600)   # bins 2000-5000 each
    check((a["top_start"], a["top_end"], a["bottom_start"], a["bottom_end"]) == (2000, 5000, 2000, 5000),
          "each side snapped to bins 2000-5000")
    check(a["top"].replace("-", "").replace(" ", "") == asm1[2000:5000], "the whole asm1 stretch on top")
    check(a["bottom"].replace("-", "").replace(" ", "") == asm2[2000:5000], "the whole asm2 stretch below")
    check(a["mismatches"] == 2 and a["ins"] == 12 and a["ins_opens"] == 1 and a["del"] == 0,
          f"exactly the planted variants: {a['mismatches']} mm, {a['ins']} ins in "
          f"{a['ins_opens']}, {a['del']} del")
    check(a["outside_top"] > 300 and a["outside_bottom"] > 300,
          f"the offset shows as outside columns: {a['outside_top']} top, {a['outside_bottom']} bottom")


def test_shift_click_in_a_repeat():
    """A degenerate rectangle aligns the single bins, and finds the true counterpart
    even 5 bp from the end of its asm1 bin inside a tandem array."""
    rng = random.Random(3)
    unit = rand_seq(16, rng)
    array = "".join(unit if rng.random() > 0.3 else unit[:5] + rand_seq(1, rng) + unit[6:]
                    for _ in range(400))
    asm1 = rand_seq(2000, rng) + array + rand_seq(2000, rng)
    off = 3300
    asm2 = rand_seq(off, rng) + asm1 + rand_seq(1000, rng)
    rd, f1, f2 = _setup(asm1, asm2, binsize=800)
    y = 4000 + 795
    a = align_rect(rd, f1, f2, y + off, y + off, y, y)
    check(a["counterpart"] == [4000 + off, 4800 + off],
          f"lands on the counterpart: {a['counterpart']}, expected {[4000 + off, 4800 + off]}")
    check(a["mismatches"] == 0 and a["ins"] == 0 and a["del"] == 0, "no differences where they overlap")


def test_reverse_orientation():
    rng = random.Random(5)
    asm1 = rand_seq(5000, rng)
    asm2 = revcomp(asm1)
    rd, f1, f2 = _setup(asm1, asm2, binsize=500)
    a = align_rect(rd, f1, f2, x0=5000 - 2250, x1=5000 - 2250, y0=2250, y1=2250)
    check(a["strand"] == "-", f"inversion detected: strand {a['strand']}")
    check(a["identity"] == 1.0 and a["compared"] == 500, f"identity {a['identity']} over {a['compared']}")
    check(a["top"] == asm1[2000:2500], "the top row stays asm1 forward")
    check(a["bottom"] == revcomp(asm2[2500:3000]), "the bottom row is the asm2 bin reverse-complemented")


def test_too_large_is_refused():
    rng = random.Random(6)
    asm1 = rand_seq(MAX_SPAN + 5000, rng)
    rd, f1, f2 = _setup(asm1, asm1, binsize=1000)
    try:
        align_rect(rd, f1, f2, 0, 10, 0, MAX_SPAN + 4000)
        check(False, "a side over MAX_SPAN must be refused")
    except ValueError as exc:
        check("at most" in str(exc), f"message says why: {exc}")


def test_long_stretch_stays_fast():
    rng = random.Random(7)
    asm1 = rand_seq(70_000, rng)
    asm2 = plant(asm1, 30_000, 12, (10_000, 50_000))
    rd, f1, f2 = _setup(asm1, asm2, binsize=1000)
    t = time.time()
    a = align_rect(rd, f1, f2, 5_000, 65_000, 5_000, 65_000)
    dt = time.time() - t
    check(a["ins"] == 12 and a["ins_opens"] == 1 and a["mismatches"] == 2,
          f"60 kb: {a['mismatches']} mm, {a['ins']} ins in {a['ins_opens']}")
    check(dt < 10, f"60 kb a side took {dt:.1f} s")


def main():
    tests = [test_layout_rows_and_counts, test_polish_joins_a_scattered_gap, test_gaps_are_left_aligned,
             test_span_at_snaps_to_bins_and_keeps_one_segment,
             test_whole_stretches_with_planted_variants, test_shift_click_in_a_repeat,
             test_reverse_orientation, test_too_large_is_refused, test_long_stretch_stays_fast]
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
