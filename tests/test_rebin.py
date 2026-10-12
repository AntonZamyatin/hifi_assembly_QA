"""The marginals' zoom-dependent binning: region.marginal_bins (Python, first frame)
and glyphs.REBIN_JS (the page, every zoom) must agree, and both must match the
server's fixed-size bins at the same bin size.

The JS half runs under node; it is skipped (and said so) when node is missing.

Run: venv/bin/python tests/test_rebin.py
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from fixtures import make_region_data

from qaview.glyphs import REBIN_JS, class_stacking
from qaview.region import MARGIN_BINS, Placements, arrange, coverage_gaps, marginal_bins, nice_bin

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def test_matches_server_bins(tmp):
    """At the server's bin size over the whole locus: the same bins, band by band."""
    rd = make_region_data(tmp, "sb")
    st = class_stacking(rd)
    for asm, cov in ((rd.asm1, rd.cov1), (rd.asm2, rd.cov2)):
        ext = rd.plot_h if asm == rd.asm1 else rd.plot_w
        b = marginal_bins(rd.segs_info(asm), rd.placements[asm], st.read_band, len(st.keys),
                          0, ext, rd.binsize)
        check(np.allclose(b["series"], np.array(st.series[asm], float)),
              f"{asm}: bins differ from the server's")
        left, _ = rd.bin_edges(asm)
        check(list(b["left"]) == list(left), f"{asm}: bin edges differ")


def test_nice_bin():
    check([nice_bin(s, 500) for s in (10, 500, 501, 1000, 1001, 2600, 7e6)] ==
          [1, 1, 2, 2, 5, 10, 20000], "1-2-5 steps, never below 1 bp")


def test_coverage_gaps():
    """Dropouts: the stretches of each segment no placement covers, exact."""
    segs = [("c", 0, 100, 0), ("d", 50, 80, 100)]          # axis 0-100, then 100-130
    p0, p1, seg = [10, 30, 25, 100], [20, 60, 40, 130], [0, 0, 0, 1]
    pl = Placements(p0=np.array(p0, float), p1=np.array(p1, float), w=np.ones(4),
                    read=np.zeros(4, np.int32), seg=np.array(seg, np.int32),
                    first=np.zeros(1, np.int32), last=np.full(1, 4, np.int32))
    left, right = coverage_gaps(segs, pl)
    check(list(zip(left, right)) == [(0, 10), (20, 25), (60, 100)],
          f"gaps {list(zip(left, right))}: before, between and after, none in a full segment")


def _js(rd, a, b, gs):
    """Run qvRebin for asm1 over [a, b] and qvBinsOf for reads `gs` under node."""
    st = class_stacking(rd)
    pl = rd.placements[rd.asm1]
    segs = rd.asm1_segs_info
    data = dict(
        iv=dict(p0=pl.p0.tolist(), p1=pl.p1.tolist(), w=pl.w.tolist(),
                read=pl.read.tolist(), seg=pl.seg.tolist()),
        segs=dict(off=[x[3] for x in segs], len=[x[2] - x[1] for x in segs],
                  contig=[x[0] for x in segs], start=[x[1] for x in segs],
                  flip=[int(f) for f in rd.flips_of(rd.asm1)]),
        reads=dict(band_class=list(st.read_band), first_y=pl.first.tolist(), last_y=pl.last.tolist()),
        n=len(st.keys), a=a, b=b, gs=gs, target=MARGIN_BINS)
    script = REBIN_JS + """
    const D = JSON.parse(require("fs").readFileSync(0, "utf8"));
    const M = {iv: {data: D.iv}, segs: {data: D.segs}, hit: {data: {}},
               areas: {class: Array.from({length: D.n}, () => ({data: {}}))},
               first: "first_y", last: "last_y"};
    const reads = {data: D.reads};
    const peak = qvRebin(M, reads, D.a, D.b, D.target);
    M.hit.data = M.next;
    const out = {peak, left: Array.from(M.next.left), seg: Array.from(M.next.seg),
                 bin_start: M.next.bin_start, bin_end: M.next.bin_end,
                 tops: [], bins: qvBinsOf(M, reads, new Set(D.gs)).sort((x, y) => x - y)};
    for (let i = 0; i < D.n; i++) out.tops.push(Array.from(M.next["cum_class_" + i]));
    process.stdout.write(JSON.stringify(out));
    """
    r = subprocess.run(["node", "-e", script], input=json.dumps(data),
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def test_js_matches_python(tmp):
    if shutil.which("node") is None:
        print("  (skipped: no node)")
        return
    rd0 = make_region_data(tmp, "js")
    n = len(rd0.asm1_segs_info)
    # as built, and reversed with every segment flipped (region.arrange)
    for rd in (rd0, arrange(rd0, "y", list(range(n))[::-1], [True] * n)):
        _js_matches(rd)


def _js_matches(rd):
    st = class_stacking(rd)
    pl = rd.placements[rd.asm1]
    for a, b in ((0, rd.plot_h), (1234.5, 3456.7), (2999, 3001), (rd.plot_h - 700, rd.plot_h)):
        got = _js(rd, a, b, [int(pl.read[0])])
        ref = marginal_bins(rd.asm1_segs_info, pl, st.read_band, len(st.keys), a, b,
                            flips=rd.flips_of(rd.asm1))
        check(got["bin_start"] == list(ref["bin_start"]) and got["bin_end"] == list(ref["bin_end"]),
              f"[{a}, {b}]: bins' genome coordinates differ")
        tops = np.cumsum(ref["series"], axis=0)
        check(got["left"] == list(ref["left"]), f"[{a}, {b}]: bin edges differ")
        check(got["seg"] == list(ref["seg"]), f"[{a}, {b}]: bin segments differ")
        check(np.allclose(np.array(got["tops"]), tops), f"[{a}, {b}]: band tops differ")
        check(abs(got["peak"] - float(tops[-1].max(initial=0))) < 1e-9, f"[{a}, {b}]: peak")
        # qvBinsOf: the bins read 0 has a placement in, by brute force.
        g = int(pl.read[0])
        want = sorted({j for i in range(pl.first[g], pl.last[g])
                       for j in range(len(ref["left"]))
                       if ref["left"][j] < pl.p1[i] and ref["right"][j] > pl.p0[i]})
        check(got["bins"] == want, f"[{a}, {b}]: bins of read {g}: {got['bins']} != {want}")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_rebin_tests_")
    tests = [(test_matches_server_bins, True), (test_nice_bin, False), (test_coverage_gaps, False),
             (test_js_matches_python, True)]
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
