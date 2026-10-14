"""Sorting one axis's segments along the other's diagonal, with flips
(region.diagonal_order / region.arrange), and everything that maps the axis
back to the genome under it.

Run: venv/bin/python tests/test_arrange.py
"""

import os
import sys
import tempfile
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from fixtures import make_region_data

from qaview.align import span_at
from qaview.glyphs import _background_index, contig_end_marks
from qaview.region import arrange, coverage_gaps, diagonal_order, marginal_bins
from qaview.remember import place

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _genome(segs, flips, p, k):
    c, s, e, off = segs[k]
    return e - (p - off) if flips[k] else s + (p - off)


def test_diagonal_order(tmp):
    """Y segments A, B, C whose reads land on X at 300-400 (backwards), 0-200 and
    200-300: laid out B, C, A, with A flipped."""
    rd = make_region_data(tmp, "do")
    segs = [("y1", 0, 100, 0), ("y2", 0, 200, 100), ("y3", 0, 100, 300)]
    # (y0, y1, x0, x1) per line, one read each
    rows = [(0, 100, 400, 300), (100, 300, 0, 200), (300, 400, 200, 300), (310, 320, 50, 60)]
    y0, y1, x0, x1 = (np.array(c, float) for c in zip(*rows))
    lines = replace(rd.lines, x0=x0, x1=x1, y0=y0, y1=y1,
                    read_idx=np.arange(4, dtype=np.int32), layer=np.zeros(4, np.int8))
    rd2 = replace(rd, asm1_segs_info=segs, lines=lines, read_names=["a", "b", "c", "d"])
    order, flips = diagonal_order(rd2, "y")
    check(order == [1, 2, 0], f"order {order}")
    check(flips == [True, False, False], f"flips {flips}")


def test_arrange_keeps_the_genome(tmp):
    """Reversed and flipped: every line end, placement and gap is the same genome
    stretch it was."""
    rd = make_region_data(tmp, "ar")
    n = len(rd.asm1_segs_info)
    order, flips = list(range(n))[::-1], [True] * n
    ra = arrange(rd, "y", order, flips)
    old, new = rd.asm1_segs_info, ra.asm1_segs_info
    fo, fn = rd.flips_of(rd.asm1), ra.flips_of(ra.asm1)
    check([s[:3] for s in new] == [old[j][:3] for j in order], "segments in the new order")
    check(new[0][3] == 0 and all(new[i + 1][3] == new[i][3] + new[i][2] - new[i][1]
                                 for i in range(n - 1)), f"contiguous offsets {new}")
    check(fn == [True] * n, f"flips in the new order {fn}")
    check(np.array_equal(ra.lines.x0, rd.lines.x0), "the other axis is untouched")

    def ends(r, segs, fl, y0, y1):
        offs = np.array([s[3] for s in segs], float)
        k = np.clip(np.searchsorted(offs, (y0 + y1) / 2, "right") - 1, 0, len(segs) - 1)
        return sorted((_genome(segs, fl, y0, k), _genome(segs, fl, y1, k)))

    for i in range(len(rd.lines)):
        a = ends(rd, old, fo, rd.lines.y0[i], rd.lines.y1[i])
        b = ends(ra, new, fn, ra.lines.y0[i], ra.lines.y1[i])
        check(np.allclose(a, b), f"line {i}: genome {a} became {b}")

    po, pn = rd.placements[rd.asm1], ra.placements[ra.asm1]
    for i in range(len(po)):
        a = sorted(_genome(old, fo, p, po.seg[i]) for p in (po.p0[i], po.p1[i]))
        b = sorted(_genome(new, fn, p, pn.seg[i]) for p in (pn.p0[i], pn.p1[i]))
        check(np.allclose(a, b) and pn.p0[i] <= pn.p1[i], f"placement {i}: {a} became {b}")
    go, gn = coverage_gaps(old, po), coverage_gaps(new, pn)
    check(abs((go[1] - go[0]).sum() - (gn[1] - gn[0]).sum()) < 1e-9, "same dropout length")

    # bins name the genome they cover, backwards on a flipped segment
    b = marginal_bins(new, pn, [0] * len(rd.read_names), 1, 0, rd.plot_h, 500, fn)
    for j in range(len(b["left"])):
        k = b["seg"][j]
        g = sorted((_genome(new, fn, b["left"][j], k), _genome(new, fn, b["right"][j], k)))
        check([b["bin_start"][j], b["bin_end"][j]] == g, f"bin {j}: {b['bin_start'][j]}-"
              f"{b['bin_end'][j]} vs {g}")
    tot_o = marginal_bins(old, po, [0] * len(rd.read_names), 1, 0, rd.plot_h, 1)["series"].sum()
    tot_n = marginal_bins(new, pn, [0] * len(rd.read_names), 1, 0, rd.plot_h, 1, fn)["series"].sum()
    check(abs(tot_o - tot_n) < 1e-6, f"coverage at 1 bp sums alike: {tot_o} vs {tot_n}")
    check(ra.genome_at(ra.asm1, new[0][3]) == (new[0][0], new[0][2] - 1),
          "a flipped segment starts at its last base")


def test_axis_consumers():
    segs = [("c", 1000, 2000, 0)]
    # the aligner: the same stretch, read either way
    fwd = span_at(segs, 100, 100, 350)
    rev = span_at(segs, 100, 650, 900, [True])
    check(fwd[3:5] == rev[3:5] == (1100, 1400), f"aligner spans {fwd[3:5]} / {rev[3:5]}")
    # remembered rectangles: placed mirrored
    meta = {"id": "m", "note": "", "created": "",
            "top": {"assembly": "A", "contig": "c", "start": 1100, "end": 1400},
            "bottom": {"assembly": "B", "contig": "c", "start": 1100, "end": 1400}}
    r = place([meta], segs, segs, "A", "B", [True], None)[0]
    check((r["y0"], r["y1"], r["x0"], r["x1"]) == (600, 900, 100, 400), f"placed {r}")


def test_end_marks_and_background(tmp):
    rd = make_region_data(tmp, "em")
    n = len(rd.asm1_segs_info)
    ra = arrange(rd, "y", list(range(n)), [True] * n)
    got = {m[:4] for m in contig_end_marks(ra, ra.asm1, tol=10 ** 9)}
    want = set()
    for c, s0, e0, off in ra.asm1_segs_info:      # flipped: contig end low, start high
        want |= {(off, 1, c, "end"), (off + e0 - s0, -1, c, "start")}
    check(got == want, f"end marks {sorted(got)} != {sorted(want)}")
    bs = 100
    nbin = -(-rd.plot_h // bs)
    check(list(_background_index(rd, rd.asm1, nbin, bs)) == list(range(nbin)),
          "the built layout maps the image onto itself")
    idx = _background_index(ra, ra.asm1, nbin, bs)
    first = ra.asm1_segs_info[0]
    k = (first[2] - first[1]) // bs
    check(list(idx[:k]) == list(range(k))[::-1], f"a flipped segment's bins run backwards: {idx[:k]}")


def test_indels_survive_sorting(tmp):
    """arrange() turns indel positions into arrays; the layout must still draw them."""
    from qaview.glyphs import report_layout
    rd = replace(make_region_data(tmp, "ind"),
                 indels={"x": [10.0, 20.0], "y": [10.0, 20.0], "read": [0, 0], "text": ["a", "b"]})
    for axis in ("x", "y"):
        n = len(rd.asm2_segs_info if axis == "x" else rd.asm1_segs_info)
        report_layout(arrange(rd, axis, list(range(n)), [False] * n))


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_arrange_tests_")
    tests = [(test_diagonal_order, True), (test_arrange_keeps_the_genome, True),
             (test_axis_consumers, False), (test_end_marks_and_background, True),
             (test_indels_survive_sorting, True)]
    for fn, needs_tmp in tests:
        print(f"{fn.__name__} ...")
        fn(tmp) if needs_tmp else fn()
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
