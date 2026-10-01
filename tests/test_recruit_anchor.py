"""Tests for the anchor extension (qatool.SRRecruiter, extension="anchor").

Synthetic data: two assemblies with the same 60 kb sequence except a block
[25,000, 35,000) that differs in asm2. Reads tile both every 1 kb, 5 kb long;
a read over the block aligns to asm2 at pid 0.995, so it favours asm1 and is not
common, while a read clear of it ties at pid 1.0 on identical sequence: common.

A seed inside the block must grow outward read by read while the outermost
crossing read is not common, and seal where it is (24,000 and 36,000), with
the matching edges of asm2 sealed by the anchors' counterparts. A coverage
gap is crossed in 1 kb steps. qaview's cache key for the legacy extension is
the one it had before `extension` existed.

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


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    pysam.faidx(path)
    return path


def make_data(tmp, reads_at, contigs1, contigs2, differs=None):
    """FASTAs, BAMs and the aln table for reads at `reads_at` [(start, length)].

    contigsN: {name: sequence}; reads are placed at the same coordinates on the
    first contig of each. `differs(start, end)` -> True if asm2 differs there.
    """
    c1, c2 = next(iter(contigs1)), next(iter(contigs2))
    fa1 = write_fasta(os.path.join(tmp, "a1.fa"), contigs1)
    fa2 = write_fasta(os.path.join(tmp, "a2.fa"), contigs2)
    bam1_reads, bam2_reads, rows = [], [], []
    for start, length in reads_at:
        name = f"r{start // 1000}"
        nm2 = 25 if differs and differs(start, start + length) else 0
        bam1_reads.append((name, c1, start, length, 0, 0))
        bam2_reads.append((name, c2, start, length, nm2, 0))
        rows.append({"assembly": "A", "read_name": name, "contig": c1,
                     "start": start, "end": start + length, "pid": 1.0})
        rows.append({"assembly": "B", "read_name": name, "contig": c2,
                     "start": start, "end": start + length, "pid": 1 - nm2 / length})
    bam1 = write_bam(os.path.join(tmp, "a1.bam"), bam1_reads,
                     {c: len(s) for c, s in contigs1.items()})
    bam2 = write_bam(os.path.join(tmp, "a2.bam"), bam2_reads,
                     {c: len(s) for c, s in contigs2.items()})
    db = write_db(os.path.join(tmp, "aln.sqlite"), rows)
    return dict(db_path=db, asm1="A", asm2="B", asm1_bam_path=bam1, asm2_bam_path=bam2,
                asm1_fasta_path=fa1, asm2_fasta_path=fa2)


def recruiter(files, extension="anchor", merge_gap=1_000):
    from qatool.SRRecruiter import Recruiter

    return Recruiter(**files, pid_threshold=0.99, plen_threshold=0.5, best_pid_delta=0.001,
                     merge_gap=merge_gap, extension_window=100, ext_max=300_000,
                     ext_rounds=100, max_rounds=100, out_dir=None, verbose=False,
                     extension=extension)


def spans(sc, asm):
    return sorted((s.contig, s.start, s.end) for s in sc if s.assembly == asm)


def test_anchor_seals_at_common_reads():
    from qatool.SegmentStorage import Segment

    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    other = list(seq)
    for i in range(BLOCK[0], BLOCK[1], 200):              # asm2 differs inside the block
        other[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[other[i]]
    differs = lambda s, e: s < BLOCK[1] and e > BLOCK[0]
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, [(i * 1000, READ) for i in range((L - READ) // 1000 + 1)],
                          {"c1": seq}, {"c2": "".join(other)}, differs)
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "c1", 28_000, 32_000)])
        check(spans(sc, "A") == [("c1", 24_000, 36_000)], f"asm1 grows to the anchors: {spans(sc, 'A')}")
        check(spans(sc, "B") == [("c2", 24_000, 36_000)], f"asm2 stops at their counterparts: {spans(sc, 'B')}")
        check(reads == {f"r{i}" for i in range(24, 32)}, f"reads inside: {sorted(reads)}")
        check(r.boundary_reads == {f"r{i}" for i in (20, 21, 22, 23, 32, 33, 34, 35)},
              f"reads crossing a sealed edge, kept apart: {sorted(r.boundary_reads)}")
        check(not (r.boundary_reads & reads), "a boundary read is never a locus read")
        log = [(e["assembly"], e["side"], e["from"], e["to"], e["why"]) for e in r.edge_log]
        check(("A", "L", 28_000, 24_000, "anchor") in log and ("A", "R", 32_000, 36_000, "anchor") in log,
              f"asm1's edges stop at anchors: {log}")
        check(r._edge_seals.get(("B", "c2", "L")) == {24_000}
              and r._edge_seals.get(("B", "c2", "R")) == {36_000},
              f"the anchors seal asm2 at the same places: {r._edge_seals}")

        # the legacy extension still runs on the same data and keeps no boundary reads
        old = recruiter(files, "legacy")
        old.run([Segment("A", "c1", 28_000, 32_000)])
        check(old.boundary_reads == set() and old.edge_log == [], "legacy is untouched")


def test_seal_inside_a_merged_segment_releases_its_reads():
    """Two loci sealed 6 kb apart merge (merge_gap 10 kb): the seals they face end up
    inside one segment, and the reads crossing them -- the reads covering the merge
    gap -- are ordinary reads of the locus again, not boundary reads."""
    from qatool.SegmentStorage import Segment

    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    other = list(seq)
    blocks = [BLOCK, (45_000, 48_000)]
    for b0, b1 in blocks:
        for i in range(b0, b1, 200):
            other[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[other[i]]
    differs = lambda s, e: any(s < b1 and e > b0 for b0, b1 in blocks)
    with tempfile.TemporaryDirectory() as tmp:
        files = make_data(tmp, [(i * 1000, READ) for i in range((L - READ) // 1000 + 1)],
                          {"c1": seq}, {"c2": "".join(other)}, differs)
        r = recruiter(files, merge_gap=10_000)
        sc, reads, _ = r.run([Segment("A", "c1", 28_000, 32_000), Segment("A", "c1", 46_000, 47_000)])
        # alone they would seal at 24,000-36,000 and 42,000-51,000
        check(spans(sc, "A") == [("c1", 24_000, 51_000)], f"merged: {spans(sc, 'A')}")
        check(reads == {f"r{i}" for i in range(24, 47)},
              f"every read inside the merged segment is recruited: missing "
              f"{sorted({f'r{i}' for i in range(24, 47)} - reads)}")
        check(r.boundary_reads == {f"r{i}" for i in (20, 21, 22, 23, 47, 48, 49, 50)},
              f"boundary reads only at the final edges: {sorted(r.boundary_reads)}")


def test_gap_is_stepped_over():
    from qatool.SegmentStorage import Segment

    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(30_000))
    with tempfile.TemporaryDirectory() as tmp:
        # reads only over the first 10 kb; the seed sits in the empty rest
        files = make_data(tmp, [(i * 1000, READ) for i in range(6)], {"g1": seq}, {"g2": seq})
        r = recruiter(files)
        sc, reads, _ = r.run([Segment("A", "g1", 20_000, 21_000)])
        check(spans(sc, "A") == [("g1", 9_000, 30_000)],
              f"left: 1 kb steps until a read crosses (r5 at 9,000), right: to the contig end: {spans(sc, 'A')}")
        log = {(e["side"], e["why"]) for e in r.edge_log}
        check(log == {("L", "anchor"), ("R", "contig_end")}, f"edge log: {r.edge_log}")
        check(reads == set() and r.boundary_reads == {"r5"}, f"{reads}, {r.boundary_reads}")


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


def main():
    tests = [test_anchor_seals_at_common_reads, test_seal_inside_a_merged_segment_releases_its_reads,
             test_gap_is_stepped_over,
             test_qaview_legacy_key_unchanged]
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
