"""Tests for the anchor extension (qatool.SRRecruiter, extension="anchor").

Synthetic data: two assemblies with the same sequence except blocks where asm2
differs. Reads tile both every 1 kb, 5 kb long, so none bridges another; a read
over a block aligns to asm2 at pid 0.995 and is not common, one clear of it ties
at pid 1.0 on identical sequence: common, and unbridged, so an anchor.

An edge must move to the nearest anchor region and take it in; the anchor's
copies -- its best loci in both assemblies -- become segments sealed on the same
outer side (the other for a reverse-complement copy) that grow inward. A read
bridging the anchor's copy in the other assembly disqualifies it. Common reads
add no segments of their own. qaview's cache key for the legacy extension is the
one it had before `extension` existed.

Run: venv/bin/python tests/test_recruit_anchor.py
"""

import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam

from fixtures import write_bam, write_db

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


L = 60_000
BLOCK = (25_000, 35_000)
READ = 5_000
SWAP = {"A": "C", "C": "G", "G": "T", "T": "A"}


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    pysam.faidx(path)
    return path


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def with_blocks(seq, blocks):
    out = list(seq)
    for b0, b1 in blocks:
        for i in range(b0, b1, 200):
            out[i] = SWAP[out[i]]
    return "".join(out)


def tile(n, blocks=(), c1="c1", c2="c2"):
    """Reads at i*1000 for i < n, the same place in both assemblies: (name, asm, contig,
    start, length, nm, flag); 25 mismatches in asm2 over a block."""
    out = []
    for i in range(n):
        s = i * 1000
        nm2 = 25 if any(s < b1 and s + READ > b0 for b0, b1 in blocks) else 0
        out += [(f"r{i}", "A", c1, s, READ, 0, 0), (f"r{i}", "B", c2, s, READ, nm2, 0)]
    return out


def _aln_type(flag, extra):
    if flag & 2048:
        return "supplementary"
    if flag & 256:
        return "secondary"
    return "primary_wsup" if "SA" in extra.get("tags", {}) else "primary"


def make_data(tmp, alns, contigs1, contigs2, extra_rows=()):
    """FASTAs, BAMs and the aln table for alignments (name, asm, contig, start, length,
    nm, flag[, extra]); flag 16 = reverse strand, 256 = secondary, 2048 =
    supplementary; `extra` as fixtures.write_bam takes it (clips for a split piece,
    an SA tag on the primary it belongs to). `extra_rows` go to the table only."""
    fa1 = write_fasta(os.path.join(tmp, "a1.fa"), contigs1)
    fa2 = write_fasta(os.path.join(tmp, "a2.fa"), contigs2)
    rows = []
    for name, asm, contig, start, length, nm, flag, *rest in alns:
        extra = rest[0] if rest else {}
        clip = sum(extra.get("clip", (0, 0)))
        rows.append({"assembly": asm, "read_name": name, "contig": contig, "start": start,
                     "end": start + length, "pid": 1 - nm / length,
                     "strand": 0 if flag & 16 else 1, "aln_type": _aln_type(flag, extra),
                     "plen": length / (length + clip), "q_len": length + clip,
                     "hard_left": extra.get("clip", (0, 0))[0],
                     "hard_right": extra.get("clip", (0, 0))[1]})
    bam = {asm: write_bam(os.path.join(tmp, f"{asm}.bam"),
                          [(n, c, s, ln, nm, f, *rest) for n, a, c, s, ln, nm, f, *rest in alns
                           if a == asm],
                          {c: len(s) for c, s in contigs.items()})
           for asm, contigs in (("A", contigs1), ("B", contigs2))}
    db = write_db(os.path.join(tmp, "aln.sqlite"), rows + list(extra_rows))
    return dict(db_path=db, asm1="A", asm2="B", asm1_bam_path=bam["A"], asm2_bam_path=bam["B"],
                asm1_fasta_path=fa1, asm2_fasta_path=fa2)


def recruiter(files, extension="anchor", merge_gap=1_000, plen=0.5):
    from qatool.SRRecruiter import Recruiter

    return Recruiter(**files, pid_threshold=0.99, plen_threshold=plen, best_pid_delta=0.001,
                     merge_gap=merge_gap, extension_window=100, ext_max=300_000,
                     ext_rounds=100, max_rounds=100, verbose=False,
                     extension=extension)


def spans(sc, asm):
    return sorted((s.contig, s.start, s.end) for s in sc if s.assembly == asm)


def test_edges_move_to_the_nearest_anchors():
    from qatool.SegmentStorage import Segment

    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    n = (L - READ) // 1000 + 1
    # r20 also aligns at 50-55 kb in both, below its best but above the thresholds:
    # a segment, as every qualifying alignment of a recruited read is
    far = [{"assembly": asm, "read_name": "r20", "contig": c, "start": 50_000, "end": 55_000,
            "pid": 0.995} for asm, c in (("A", "c1"), ("B", "c2"))]
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, tile(n, [BLOCK]), {"c1": seq}, {"c2": with_blocks(seq, [BLOCK])}, far)
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "c1", 28_000, 32_000)])
        # r21..r34 overlap the block; r20 (20-25 kb) and r35 (35-40 kb) are the anchors
        check(spans(sc, "A")[0] == ("c1", 20_000, 40_000), f"asm1 takes in its anchors: {spans(sc, 'A')}")
        check(spans(sc, "B")[0] == ("c2", 20_000, 40_000), f"asm2 from the anchors' copies: {spans(sc, 'B')}")
        for asm in ("A", "B"):
            check(any(s <= 50_000 and e >= 55_000 for _, s, e in spans(sc, asm)[1:]),
                  f"{asm}: r20's other qualifying alignment is a segment too: {spans(sc, asm)}")
        check(reads >= {f"r{i}" for i in range(20, 36)}, f"reads inside, anchors too: {sorted(reads)}")
        log = [(e["assembly"], e["side"], e["from"], e["to"], e["why"], e["read"]) for e in r.edge_log]
        check(("A", "L", 28_000, 20_000, "anchor", "r20") in log
              and ("A", "R", 32_000, 40_000, "anchor", "r35") in log, f"asm1's edges: {log}")
        check(("B", "L", None, 20_000, "copy", "r20") in log
              and ("B", "R", None, 40_000, "copy", "r35") in log, f"their copies in asm2: {log}")
        check(20_000 in r._edge_seals.get(("B", "c2", "L"), ())
              and 40_000 in r._edge_seals.get(("B", "c2", "R"), ()), f"asm2's seals: {r._edge_seals}")

        old = recruiter(files, "legacy")
        sc, _, _ = old.run([Segment("A", "c1", 28_000, 32_000)])
        check(old.edge_log == [] and ("c1", 50_000, 55_000) in spans(sc, "A"),
              f"legacy is untouched, and still follows r20 to 50 kb: {spans(sc, 'A')}")


def test_anchor_bridged_in_the_other_assembly_is_not_one():
    """x aligns to asm2 only, over 19.5-25.5 kb: it bridges r20's copy there, so r20
    is no anchor although nothing bridges it in asm1; r19 (19-24 kb) is."""
    from qatool.SegmentStorage import Segment

    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    n = (L - READ) // 1000 + 1
    with tempfile.TemporaryDirectory() as tmp:
        alns = tile(n, [BLOCK]) + [("x", "B", "c2", 19_500, 6_000, 0, 0)]
        files = make_data(tmp, alns, {"c1": seq}, {"c2": with_blocks(seq, [BLOCK])})
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "c1", 28_000, 32_000)])
        check(spans(sc, "A") == [("c1", 19_000, 40_000)], f"asm1 passes r20 for r19: {spans(sc, 'A')}")
        check(spans(sc, "B") == [("c2", 19_000, 40_000)], f"and asm2 with it: {spans(sc, 'B')}")
        check("x" in reads, "x is a locus read of asm2")


def test_repeat_copies_get_segments():
    """A 6 kb repeat at 10-16 kb, reverse-complemented at 50-56 kb, and again at
    80-86 kb, in both assemblies; asm2 differs at 16-26 kb, next to the first copy.
    A read inside a copy aligns to all three. Left of the seed the nearest anchor is
    r11 (11-16 kb), in the repeat: each copy is sealed on its outer side -- the
    reverse-complement one on its right -- and grows inward to the next anchor."""
    from qatool.SegmentStorage import Segment

    rng = random.Random(5)
    size, rep = 100_000, "".join(rng.choice("ACGT") for _ in range(6_000))
    copies = [(10_000, 1), (50_000, -1), (80_000, 1)]          # (start, orientation)
    seq = list("".join(rng.choice("ACGT") for _ in range(size)))
    for c0, o in copies:
        seq[c0:c0 + 6_000] = rep if o == 1 else revcomp(rep)
    seq = "".join(seq)
    block = (16_000, 26_000)
    alns = []
    for i in range(size // 1000 - 4):
        s = i * 1000
        nm2 = 25 if s < block[1] and s + READ > block[0] else 0
        home = next(((c0, o) for c0, o in copies if c0 <= s and s + READ <= c0 + 6_000), None)
        if home is None:
            alns += [(f"r{i}", "A", "c1", s, READ, 0, 0), (f"r{i}", "B", "c2", s, READ, nm2, 0)]
            continue
        c0, o = home
        u = s - c0 if o == 1 else c0 + 6_000 - s - READ          # offset in the repeat
        for c1, o1 in copies:
            at = c1 + u if o1 == 1 else c1 + 6_000 - u - READ
            flag = (0 if o1 == o else 16) | (0 if c1 == c0 else 256)
            alns += [(f"r{i}", "A", "c1", at, READ, 0, flag), (f"r{i}", "B", "c2", at, READ, 0, flag)]
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, alns, {"c1": seq}, {"c2": with_blocks(seq, [block])})
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "c1", 20_000, 22_000)])
        want = [("c1", 11_000, 31_000), ("c1", 49_000, 55_000), ("c1", 81_000, 87_000)]
        check(spans(sc, "A") == want, f"the locus and the repeat's other copies in asm1: {spans(sc, 'A')}")
        check(spans(sc, "B") == [(("c2",) + w[1:]) for w in want], f"and in asm2: {spans(sc, 'B')}")
        check(r._edge_seals.get(("A", "c1", "R"), set()) >= {55_000}
              and 50_000 not in r._edge_seals.get(("A", "c1", "L"), set()),
              f"the reverse-complement copy is sealed on its right: {r._edge_seals}")
        check({"r11", "r49", "r50", "r81", "r82"} <= reads, f"anchors of the copies: {sorted(reads)}")


def test_merged_loci_recruit_every_read_inside():
    """Two loci whose anchors meet at 40 kb merge; every read inside is recruited."""
    from qatool.SegmentStorage import Segment

    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    blocks = [BLOCK, (45_000, 48_000)]
    n = (L - READ) // 1000 + 1
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, tile(n, blocks), {"c1": seq}, {"c2": with_blocks(seq, blocks)})
        r = recruiter(files, merge_gap=10_000)
        sc, reads, _ = r.run([Segment("A", "c1", 28_000, 32_000), Segment("A", "c1", 46_000, 47_000)])
        check(spans(sc, "A") == [("c1", 20_000, 53_000)], f"merged: {spans(sc, 'A')}")
        check(reads == {f"r{i}" for i in range(20, 49)},
              f"every read inside: {sorted(set(f'r{i}' for i in range(20, 49)) ^ reads)}")


def test_gap_is_crossed_to_the_nearest_anchor():
    from qatool.SegmentStorage import Segment

    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(30_000))
    with tempfile.TemporaryDirectory() as tmp:
        # reads only over the first 10 kb; the seed sits in the empty rest
        files = make_data(tmp, tile(6, c1="g1", c2="g2"), {"g1": seq}, {"g2": seq})
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "g1", 20_000, 21_000)])
        check(spans(sc, "A") == [("g1", 5_000, 30_000)],
              f"left: to r5 (5-10 kb) across the gap, right: to the contig end: {spans(sc, 'A')}")
        check(spans(sc, "B") == [("g2", 5_000, 30_000)], f"asm2 from r5's copy: {spans(sc, 'B')}")
        log = {(e["assembly"], e["side"], e["why"]) for e in r.edge_log}
        check(log == {("A", "L", "anchor"), ("A", "R", "contig_end"), ("B", "L", "copy"),
                      ("B", "R", "contig_end")}, f"edge log: {r.edge_log}")
        check(reads == {"r5"}, f"{reads}")


def test_qaview_legacy_key_unchanged():
    import hashlib
    import json

    from qaview.pipeline import RECRUIT_VERSION, recruit_key
    from qaview.project import RECRUITMENT_DEFAULTS, Project, Seed, validate_recruitment

    seed = Seed("s", "bed", "A", "c1", 10, 20)
    p = Project(path="p.toml", name="p", asm1=type("C", (), {"name": "A", "fasta": "a", "bam": "b"})(),
                asm2=type("C", (), {"name": "B", "fasta": "c", "bam": "d"})(), db="db",
                library_assembly="B", seeds=[seed], recruitment=dict(RECRUITMENT_DEFAULTS),
                resources={}, output_dir="out", seed_sources=[])
    before = {k: v for k, v in RECRUITMENT_DEFAULTS.items() if k != "extension"}
    old = hashlib.sha1(json.dumps({
        "version": RECRUIT_VERSION, "db": "db", "asm": [["A", "a", "b"], ["B", "c", "d"]],
        "seed": ["A", "c1", 10, 20], "params": before}, sort_keys=True).encode()).hexdigest()[:16]
    check(recruit_key(p, seed) == old, "a legacy recruitment keeps its cache key")
    check(recruit_key(p, seed, {**p.recruitment, "extension": "anchor"}) != old,
          "an anchor recruitment has its own")
    check(validate_recruitment({"extension": "anchor"}) == {"extension": "anchor"}, "anchor is valid")
    try:
        validate_recruitment({"extension": "nope"})
        check(False, "an unknown extension is refused")
    except ValueError as exc:
        check("one of legacy, anchor" in str(exc), str(exc))


def test_common_read_brings_its_other_copies():
    """x is common: identical at 20.5-24.5 kb, inside the locus, and at 50-54 kb, far
    outside it, in both assemblies. Its other copy must become a segment in both,
    or the scorer would see it in one candidate only (+/-C)."""
    from qatool.SegmentStorage import Segment

    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    seq = seq[:50_000] + seq[20_500:24_500] + seq[54_000:]      # the copy
    n = (L - READ) // 1000 + 1
    alns = tile(n, [BLOCK])
    for asm, c in (("A", "c1"), ("B", "c2")):
        alns += [("x", asm, c, 20_500, 4_000, 0, 0), ("x", asm, c, 50_000, 4_000, 0, 256)]
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, alns, {"c1": seq}, {"c2": with_blocks(seq, [BLOCK])})
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "c1", 28_000, 32_000)])
        for asm in ("A", "B"):
            sp = spans(sc, asm)
            check(any(s <= 50_000 and e >= 54_000 for _, s, e in sp),
                  f"{asm}: x's copy at 50-54 kb is a segment: {sp}")
        check("x" in reads, "x is a locus read")


def test_split_alignments_take_no_part_in_recruitment():
    """Reads split by a supplementary alignment, whole alignments at plen 0.999:
    s is whole in A inside the locus, split in B into 26-28.5 kb and 52-54.5 kb;
    t is split in A itself, a piece inside the locus (33-35.5 kb) and one at
    56-58.5 kb. Partial alignments recruit nothing -- random ones would tie in
    unrelated regions: s is recruited whole and brings no piece, t not at all."""
    from qatool.SegmentStorage import Segment

    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    n = (L - READ) // 1000 + 1
    half = {"clip": (0, 2_500)}, {"clip": (2_500, 0)}
    sa = {"tags": {"SA": "x,1,+,1M,60,0;"}}
    alns = tile(n, [BLOCK]) + [
        ("s", "A", "c1", 26_000, 5_000, 0, 0),
        ("s", "B", "c2", 26_000, 2_500, 0, 0, {**half[0], **sa}),
        ("s", "B", "c2", 52_000, 2_500, 0, 2048, half[1]),
        ("t", "A", "c1", 33_000, 2_500, 0, 0, {**half[0], **sa}),
        ("t", "A", "c1", 56_000, 2_500, 0, 2048, half[1]),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, alns, {"c1": seq}, {"c2": with_blocks(seq, [BLOCK])})
        sc, reads, _ = recruiter(files, plen=0.999).run([Segment("A", "c1", 28_000, 32_000)])
    covers = lambda asm, lo, hi: any(s < hi and e > lo for _, s, e in spans(sc, asm))
    check("s" in reads and "t" not in reads, f"s whole, t only split: {sorted(reads & {'s', 't'})}")
    check(not covers("B", 52_000, 54_500) and not covers("A", 56_000, 58_500),
          f"no piece brings a segment: {spans(sc, 'A')} {spans(sc, 'B')}")


def test_merge_gap_either_order():
    """Segments merge_gap apart merge whichever is added first."""
    from qatool.SegmentStorage import Segment, SegmentCollection
    for gap in (0, 10):
        a, b = Segment("A", "c", 0, 100), Segment("A", "c", 100 + gap, 200)
        n = [len(SegmentCollection(order, merge_gap=gap)) for order in ([a, b], [b, a])]
        check(n == [1, 1], f"gap {gap}: {n}")


def main():
    tests = [test_edges_move_to_the_nearest_anchors, test_anchor_bridged_in_the_other_assembly_is_not_one,
             test_repeat_copies_get_segments, test_merged_loci_recruit_every_read_inside,
             test_gap_is_crossed_to_the_nearest_anchor, test_common_read_brings_its_other_copies,
             test_split_alignments_take_no_part_in_recruitment,
             test_qaview_legacy_key_unchanged, test_merge_gap_either_order]
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
