"""Tests for drawing every read at its best placements (qaview.region): the dotplot
runs between a read's best placements only, whatever its class, and each marginal
splits the read's weight of 1 over its best placements in that assembly, counted
genome-wide -- a best copy outside the locus still takes its share.

Synthetic: one 20 kb contig per assembly, the locus its first 10 kb, bins of 500 bp.

Run: venv/bin/python tests/test_best_placements.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from fixtures import write_bam, write_db

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


# (read, assembly, start, end, pid, aln_type)
ALNS = [
    ("u", "A", 1000, 2000, 1.0, "primary"), ("u", "B", 1000, 2000, 1.0, "primary"),
    # two best copies in A, and a third 0.5% worse: qualifying, not best
    ("rep", "A", 3000, 4000, 1.0, "primary"), ("rep", "A", 5000, 6000, 1.0, "secondary"),
    ("rep", "A", 7000, 8000, 0.995, "secondary"), ("rep", "B", 3000, 4000, 1.0, "primary"),
    # a second best copy in A outside the locus (10-20 kb)
    ("out", "A", 8500, 9500, 1.0, "primary"), ("out", "A", 15000, 16000, 1.0, "secondary"),
    ("out", "B", 8500, 9500, 1.0, "primary"),
    # a supplementary piece of a best primary: it counts at its primary's share
    ("sup", "A", 2000, 2500, 1.0, "primary_wsup"), ("sup", "A", 6000, 6500, 0.995, "supplementary"),
    ("sup", "B", 2000, 2500, 1.0, "primary"),
]
FLAG = {"primary": 0, "primary_wsup": 0, "secondary": 256, "supplementary": 2048}


def build(tmp):
    from qatool.SegmentStorage import Segment
    from qaview.region import RegionData

    bams = {}
    for asm, contig in (("A", "cA"), ("B", "cB")):
        recs = [(r, contig, s, e - s, round((1 - pid) * (e - s)), FLAG[t])
                for r, a, s, e, pid, t in ALNS if a == asm]
        bams[asm] = write_bam(os.path.join(tmp, f"{asm}.bam"), recs, {contig: 20_000})
    rows = [{"assembly": a, "read_name": r, "contig": "cA" if a == "A" else "cB",
             "start": s, "end": e, "pid": pid, "aln_type": t} for r, a, s, e, pid, t in ALNS]
    db = write_db(os.path.join(tmp, "aln.sqlite"), rows)
    segs = [Segment("A", "cA", 0, 10_000), Segment("B", "cB", 0, 10_000)]
    return RegionData.build(segs, {"u", "rep", "out", "sup"}, db, bams["A"], bams["B"],
                            asm1="A", asm2="B", pid_threshold=0.99, plen_threshold=0.5,
                            binsize=500, best_pid_delta=0.001, pid_tie_eps=0.0001)


def test_lines_run_between_best_placements_only():
    with tempfile.TemporaryDirectory() as tmp:
        rd = build(tmp)
    names = [rd.read_names[i] for i in rd.lines.read_idx]
    count = {n: names.count(n) for n in set(names)}
    check(count.get("u") == 1, f"a unique read: one line, {count.get('u')}")
    check(count.get("rep") == 2, f"two best copies against one: two lines, not three: {count.get('rep')}")
    ys = sorted(round(min(y0, y1)) for n, y0, y1 in zip(names, rd.lines.y0, rd.lines.y1) if n == "rep")
    check(ys == [3000, 5000], f"at the two best copies, not the worse one at 7 kb: {ys}")
    check(count.get("out") == 1, f"its copy outside the locus draws nothing: {count.get('out')}")


def test_marginals_split_each_read_over_its_best_placements():
    with tempfile.TemporaryDirectory() as tmp:
        rd = build(tmp)
    a = rd.cov1.total()                      # asm A, bins of 500 bp over 0-10 kb
    at = lambda pos: a[pos // 500]
    check(np.isclose(at(1500), 1.0), f"unique read: weight 1 ({at(1500)})")
    check(np.isclose(at(3500), 0.5) and np.isclose(at(5500), 0.5),
          f"two best copies: 1/2 each ({at(3500)}, {at(5500)})")
    check(np.isclose(at(7500), 0.0), f"the worse copy adds nothing ({at(7500)})")
    check(np.isclose(at(9000), 0.5),
          f"a best copy outside the locus still takes its share: 1/2 here ({at(9000)})")
    check(np.isclose(at(2250), 1.0) and np.isclose(at(6250), 1.0),
          f"the primary and its supplementary piece both at the primary's share ({at(2250)}, {at(6250)})")
    b = rd.cov2.total()
    check(np.isclose(b.sum(), 2 * 1.0 + 2 * 1.0 + 2 * 1.0 + 1.0),
          f"asm B: each read once over its bins ({b.sum()})")


def test_one_sided_reads_show_partial_alignments_on_selection():
    """bp is whole in A only (1-6 kb) -- asm1_only -- and split in B: a primary piece
    at 1-3.5 kb carrying an SA tag and a supplementary at 7-9.5 kb (plen 0.5 each). It
    stays asm1_only and B's marginal does not count it; it gains two partial rows,
    one per piece, never drawn but when selected. cl, whole in A and only clipped in
    B at cB's end, gains one, labelled as at a contig end. d is whole in both (a
    full-length secondary in B) and split in B as well: not one-sided, no partial
    rows. Without a partial plen there are none."""
    from qatool.SegmentStorage import Segment
    from qaview.region import RegionData

    recs = {"A": [("bp", "cA", 1000, 5000, 0, 0), ("w", "cA", 6000, 3000, 0, 0),
                  ("d", "cA", 12000, 5000, 0, 0), ("cl", "cA", 14000, 4000, 0, 0),
                  ("sv", "cA", 2000, 5000, 0, 0), ("jx", "cA", 3000, 5000, 0, 0)],
            "B": [("bp", "cB", 1000, 2500, 0, 0, {"clip": (0, 2500), "tags": {"SA": "cB,7001,+,2500M,60,0;"}}),
                  ("bp", "cB", 7000, 2500, 0, 2048, {"clip": (2500, 0)}),
                  ("w", "cB", 6000, 3000, 0, 0),
                  ("d", "cB", 2000, 2500, 0, 0, {"clip": (0, 2500), "tags": {"SA": "cB,8001,+,2500M,60,0;"}}),
                  ("d", "cB", 8000, 2500, 0, 2048, {"clip": (2500, 0)}),
                  ("d", "cB", 12000, 5000, 0, 256),
                  ("cl", "cB", 17500, 2500, 0, 0, {"clip": (0, 1500)})]}
    rows = [{"assembly": "A", "read_name": "bp", "contig": "cA", "start": 1000, "end": 6000},
            {"assembly": "A", "read_name": "w", "contig": "cA", "start": 6000, "end": 9000},
            {"assembly": "B", "read_name": "bp", "contig": "cB", "start": 1000, "end": 3500,
             "aln_type": "primary_wsup", "plen": 0.5},
            {"assembly": "B", "read_name": "bp", "contig": "cB", "start": 7000, "end": 9500,
             "aln_type": "supplementary", "plen": 0.5},
            {"assembly": "B", "read_name": "w", "contig": "cB", "start": 6000, "end": 9000},
            {"assembly": "A", "read_name": "d", "contig": "cA", "start": 12000, "end": 17000},
            {"assembly": "B", "read_name": "d", "contig": "cB", "start": 2000, "end": 4500,
             "aln_type": "primary_wsup", "plen": 0.5},
            {"assembly": "B", "read_name": "d", "contig": "cB", "start": 8000, "end": 10500,
             "aln_type": "supplementary", "plen": 0.5},
            {"assembly": "B", "read_name": "d", "contig": "cB", "start": 12000, "end": 17000,
             "aln_type": "secondary"},
            {"assembly": "A", "read_name": "cl", "contig": "cA", "start": 14000, "end": 18000},
            {"assembly": "B", "read_name": "cl", "contig": "cB", "start": 17500, "end": 20000,
             "plen": 0.625, "q_len": 4000, "hard_right": 1500},
            # sv: whole-length in B but short of the display pid -- a 600 bp deletion inside
            {"assembly": "A", "read_name": "sv", "contig": "cA", "start": 2000, "end": 7000},
            {"assembly": "B", "read_name": "sv", "contig": "cB", "start": 1000, "end": 6600,
             "pid": 0.993, "cigar": "2000M600D3000M"},
            # jx: clipped in B at pid 0.993, a junction overrun -- above partial pid 0.99 only
            {"assembly": "A", "read_name": "jx", "contig": "cA", "start": 3000, "end": 8000},
            {"assembly": "B", "read_name": "jx", "contig": "cB", "start": 9000, "end": 11500,
             "pid": 0.993, "plen": 0.5, "q_len": 5000, "hard_right": 2500}]
    with tempfile.TemporaryDirectory() as tmp:
        bams = {a: write_bam(os.path.join(tmp, f"{a}.bam"), recs[a], {c: 20_000})
                for a, c in (("A", "cA"), ("B", "cB"))}
        db = write_db(os.path.join(tmp, "aln.sqlite"), rows)
        segs = [Segment("A", "cA", 0, 20_000), Segment("B", "cB", 0, 20_000)]
        kw = dict(asm1="A", asm2="B", pid_threshold=0.997, plen_threshold=0.999, binsize=500,
                  best_pid_delta=0.001, pid_tie_eps=0.0001)
        reads = {"bp", "w", "d", "cl", "sv", "jx"}
        rd = RegionData.build(segs, reads, db, bams["A"], bams["B"], supp_plen_threshold=0.3, **kw)
        relaxed = RegionData.build(segs, reads, db, bams["A"], bams["B"], supp_plen_threshold=0.3,
                                   partial_pid_threshold=0.99, min_indel=50, **kw)
        full = RegionData.build(segs, reads, db, bams["A"], bams["B"], **kw)
    names = [rd.read_names[i] for i in rd.lines.read_idx]
    part = [i for i, lay in enumerate(rd.lines.layer) if lay == 4]
    pnames = sorted(names[i] for i in part)
    check(rd.read_class.get("bp") == "asm1_only" and "breakpoint" not in rd.bucket_keys,
          f"no class of its own: {rd.read_class.get('bp')}, {rd.bucket_keys}")
    check(pnames == ["bp", "bp", "cl"] and rd.n_partial == 3, f"partial rows: {pnames}")
    check(sum(1 for i, n in enumerate(names) if n == "d" and rd.lines.layer[i] == 4) == 0,
          "a read whole in both has none")
    labels = sorted(rd.lines.part[i] for i in part)
    check(labels == ["partial primary, clipped in B, pid 1.0000, plen 0.62, at a contig end",
                     "partial primary, split in B, pid 1.0000, plen 0.50, inside a contig",
                     "partial supplementary in B, pid 1.0000, plen 0.50, inside a contig"], f"labels: {labels}")
    check(not rd.indels, "at the display pid, sv's and jx's alignments are left out")

    # at partial pid 0.99: jx's junction piece, and sv's whole-length alignment with its deletion
    rpart = [i for i, lay in enumerate(relaxed.lines.layer) if lay == 4]
    rnames = sorted(relaxed.read_names[relaxed.lines.read_idx[i]] for i in rpart)
    check(rnames == ["bp", "bp", "cl", "jx", "sv"], f"partial pid 0.99 adds jx and sv: {rnames}")
    sv = [relaxed.lines.part[i] for i in rpart if relaxed.read_names[relaxed.lines.read_idx[i]] == "sv"]
    check(sv == ["whole-length in B, pid 0.9930, plen 1.00; 600 bp deletion at cB:3,301"], f"sv: {sv}")
    ind = relaxed.indels
    check(len(ind.get("x", [])) == 1 and relaxed.read_names[ind["read"][0]] == "sv"
          and abs(ind["x"][0] - 3300) < 1 and ind["text"][0] == "600 bp deletion in B at cB:3,301",
          f"one cross, at the deletion's middle on B's axis: {ind}")
    b = rd.cov2.total()
    check(np.isclose(b[2000 // 500], 0.0) and np.isclose(b[8000 // 500], 1.0),
          f"B's marginal does not count the pieces (w alone at 8 kb): {b[4]}, {b[16]}")
    check(full.n_partial == 0, "without a partial plen there are none")

    from qaview.glyphs import _line_sources
    src, _, tables = _line_sources(rd, "class")
    d = src.data
    shown = [(float(d["alpha"][i]), float(d["sel_alpha"][i])) for i in part]
    check(all(a == 0.0 and s > 0.5 for a, s in shown), f"hidden, shown on selection: {shown}")
    check(len(tables["parts"]) == 4 and tables["parts"][0] == "full-length", f"tooltip: {tables['parts']}")
    hidden, _, _ = _line_sources(rd, "class", visible_layers={0, 1, 2, 3})
    check(all(float(hidden.data["sel_alpha"][i]) == 0.0 for i in part), "the partial layer can be switched off")


def main():
    tests = [test_lines_run_between_best_placements_only,
             test_marginals_split_each_read_over_its_best_placements,
             test_one_sided_reads_show_partial_alignments_on_selection]
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
