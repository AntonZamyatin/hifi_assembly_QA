"""Tests for the k-mer similarity background (qaview.similarity).

Synthetic: asm1 = pad + A + R + B + C; asm2 = pad + A + R + B + R1 + R2 + C, where
R1 is R with 0.2% substitutions (identity 0.998) and R2 the reverse complement of
R with 1% (0.990). Segments start past the pad, and asm2's extra copies shift its
bin grid against asm1's, so the diagonal's bins do not line up. Expected: an
unbroken diagonal near 100%; R's bins matching R1 at about 0.998 -- above the
recruitment's 0.997 -- and R2 at about 0.99, below it.

Run: venv/bin/python tests/test_similarity.py
"""

import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pysam

from qatool.SegmentStorage import Segment

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(s, rate, rng):
    out = list(s)
    for i in rng.sample(range(len(s)), int(round(rate * len(s)))):
        out[i] = rng.choice([c for c in "ACGT" if c != out[i]])
    return "".join(out)


def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n{seq}\n")
    pysam.faidx(path)
    return path


def test_canonical_kmers_match_their_reverse_complement():
    from qaview.similarity import canonical_kmers

    rng = random.Random(1)
    s = "".join(rng.choice("ACGT") for _ in range(500))
    p1, c1 = canonical_kmers(s, 21)
    p2, c2 = canonical_kmers(revcomp(s), 21)
    check(sorted(c1.tolist()) == sorted(c2.tolist()), "a sequence and its reverse complement share every canonical k-mer")
    p3, _ = canonical_kmers(s[:100] + "N" + s[101:], 21)
    check(len(p3) == len(p1) - 21, f"an N removes the 21 k-mers over it: {len(p1)} -> {len(p3)}")


def test_diagonal_and_repeat_copies(max_pairs=None):
    from qaview.similarity import MAX_PAIRS, compute

    rng = random.Random(9)
    R_ = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    pad, A, R, B, C = R_(1000), R_(30_000), R_(10_000), R_(30_000), R_(30_000)
    R1, R2 = mutate(R, 0.002, rng), revcomp(mutate(R, 0.010, rng))
    s1 = pad + A + R + B + C
    s2 = pad + A + R + B + R1 + R2 + C
    with tempfile.TemporaryDirectory() as tmp:
        f1 = write_fasta(os.path.join(tmp, "a1.fa"), "c1", s1)
        f2 = write_fasta(os.path.join(tmp, "a2.fa"), "c2", s2)
        segs = [Segment("A", "c1", 1000, len(s1)), Segment("B", "c2", 1000, len(s2))]
        res = compute(segs, "A", "B", f1, f2, binsize=2000, max_pairs=max_pairs or MAX_PAIRS)
    if max_pairs:
        check(res["settings"]["sampling_scale"] > 1, f"sampling is on: {res['settings']}")
    x0, y0, ident = res["x0"], res["y0"], res["identity"]

    # the diagonal: every full asm1 bin has a cell at its own position (A, R, B) or
    # 20 kb further in asm2 (C, after R1 and R2), at near 100%
    full_bins = list(range(0, 100_000 - 2000 + 1, 2000))
    missing = []
    for y in full_bins:
        shift = 0 if y < 70_000 else 20_000
        near = (y0 == y) & (np.abs(x0 - (y + shift)) <= 2000) & (ident >= 0.999)
        if not near.any():
            missing.append(y)
    check(not missing, f"the diagonal is unbroken despite bins out of line: missing at {missing[:10]}")

    # R (asm1 30-40 kb) against R1 (asm2 70-80 kb) and R2 (80-90 kb)
    in_R = (y0 >= 30_000) & (y0 < 40_000)
    r1 = in_R & (x0 >= 69_000) & (x0 < 81_000)
    r2 = in_R & (x0 >= 79_000) & (x0 < 91_000)
    check(r1.sum() >= 4, f"R's bins find the 99.8% copy: {r1.sum()} cells")
    check(r1.any() and np.median(ident[r1]) >= 0.997 and np.median(ident[r1]) < 0.9995,
          f"estimated at about 0.998 -- above 0.997: {np.round(ident[r1], 4)}")
    check(r2.sum() >= 4, f"and the reverse-complement copy: {r2.sum()} cells")
    check(r2.any() and 0.985 <= np.median(ident[r2]) < 0.997,
          f"estimated at about 0.990 -- below 0.997, so not drawn at the default: {np.round(ident[r2], 4)}")
    off = ~in_R & (np.abs(x0 - y0) > 4000) & ~((y0 >= 70_000) & (np.abs(x0 - y0 - 20_000) <= 4000))
    check(not off.any(), f"nothing else off the diagonal: {list(zip(x0[off][:5], y0[off][:5]))}")


def test_hash_sampling_keeps_the_picture():
    """The same with k-mers sampled 1 in ~3: estimates stay unbiased enough to keep
    the diagonal whole and the two copies either side of 0.997."""
    test_diagonal_and_repeat_copies(max_pairs=80_000)


def test_tandem_array_is_not_read_as_identical():
    """A tandem array of 1 kb units each 1% from a common unit (so ~2% between
    units): bins two units long share most k-mers with their neighbours too. Each
    k-mer must count once against a two-bin window, or every off-diagonal cell
    reads as 100%."""
    from qaview.similarity import compute

    rng = random.Random(4)
    unit = "".join(rng.choice("ACGT") for _ in range(1000))
    array = "".join(mutate(unit, 0.01, rng) for _ in range(40))
    flank = "".join(rng.choice("ACGT") for _ in range(4000))
    seq = flank + array + flank
    with tempfile.TemporaryDirectory() as tmp:
        f1 = write_fasta(os.path.join(tmp, "a1.fa"), "c1", seq)
        f2 = write_fasta(os.path.join(tmp, "a2.fa"), "c2", seq)
        res = compute([Segment("A", "c1", 0, len(seq)), Segment("B", "c2", 0, len(seq))],
                      "A", "B", f1, f2, binsize=2000)
    a, b, ident = res["a"], res["b"], res["identity"]
    inside = (a >= 3) & (a <= 20) & (b >= 3) & (b <= 20)
    off = inside & (np.abs(a - b) >= 2)
    diag = inside & (a == b)
    check(off.any() and ident[off].max() < 0.995,
          f"units ~2% apart are estimated below 0.995 off the diagonal: max {ident[off].max() if off.any() else None}")
    check(diag.sum() >= 15 and ident[diag].min() >= 0.999, f"the diagonal stays near 100%: {np.round(ident[diag], 4)}")


def test_image_shows_cells_from_the_floor_up():
    from qaview.similarity import image_rgba

    res = {"n1": 2, "n2": 3, "a": np.array([0, 1]), "b": np.array([2, 0]),
           "identity": np.array([0.999, 0.99], dtype=np.float32)}
    img = image_rgba(res, 0.997)
    check(img.shape == (2, 3) and img[0, 2] >> 24 == 255 and img[1, 0] == 0,
          f"only the cell at or above the floor is drawn: {img}")


def test_uniform_repeat_bin_does_not_light_a_row():
    """A bin of a uniform repeat (few distinct k-mers) is contained in every bin that
    also holds that unit, but they are not contained in it: taken both ways, those
    cells stay dark, while the flanks' diagonal stays (ptg000114l_translocation_001,
    where satellite bins lit whole rows)."""
    from qaview.similarity import compute

    rng = random.Random(8)
    R_ = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    # 600 bp: a bin of it keeps 30% of a full bin's distinct k-mers, above MIN_KMERS
    U, A, B = R_(600), R_(4_000), R_(4_000)
    s1 = A + U * 10 + B                                    # 6 kb of the bare unit
    s2 = A + "".join(U + R_(1_400) for _ in range(5)) + B  # the unit beside sequence of its own
    with tempfile.TemporaryDirectory() as tmp:
        f1 = write_fasta(os.path.join(tmp, "a1.fa"), "c1", s1)
        f2 = write_fasta(os.path.join(tmp, "a2.fa"), "c2", s2)
        res = compute([Segment("A", "c1", 0, len(s1)), Segment("B", "c2", 0, len(s2))],
                      "A", "B", f1, f2, binsize=2000)
    a, b = res["a"], res["b"]
    rows_u = (a >= 2) & (a <= 4)                           # asm1 4-10 kb: the bare unit
    cols_u = (b >= 2) & (b <= 6)                           # asm2 4-14 kb: unit + own sequence
    check(not (rows_u & cols_u).any(),
          f"the bare unit's bins are not drawn against them: {list(zip(a[rows_u & cols_u], b[rows_u & cols_u]))}")
    flank = ((a <= 1) & (b <= 1)) | ((a >= 5) & (b >= 7))
    check(flank.sum() >= 4 and res["identity"][flank].min() >= 0.999,
          f"the flanks' diagonal stays: {list(zip(a[flank], b[flank]))}")


def test_view_draws_from_the_recruitment_pid():
    """The background's floor is the locus's recruitment pid, not the display pid:
    lowering the display pid to see more reads must not flood the background.
    On the smallest HG002 locus, with a stand-in similarity result; skipped
    without the data."""
    import pickle
    import tempfile

    from fastapi.testclient import TestClient

    path = os.path.join(os.path.dirname(__file__), "..", "examples", "projects", "hg002.toml")
    if not os.path.exists(path):
        print("  SKIP: no examples/projects/hg002.toml")
        return
    from qaview.pipeline import Pipeline
    from qaview.project import load_project
    from qaview.server import create_app

    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_CONFIG_HOME"] = tmp
        p = load_project(path)
        pipe = Pipeline(p)
        pipe._recruit_all()
        if not pipe.loci:
            print("  SKIP: no HG002 recruitments on disk")
            return
        loc = min(pipe.loci.values(), key=lambda l: l.n_reads)
        rec_pid = loc.params["pid_threshold"]
        fake = {"n1": 2, "n2": 2, "binsize": 1000, "a": np.array([0, 1, 1]), "b": np.array([0, 1, 0]),
                "identity": np.array([rec_pid + 0.0005, rec_pid - 0.0005, 0.96], dtype=np.float32),
                "settings": {"sampling_scale": 1}}
        pipe.request_background = lambda lid: None
        pipe.background = lambda lid: fake
        with TestClient(create_app(p, pipe, start=False)) as c:
            for display in (0.99, rec_pid, 0.999):
                bg = c.get(f"/api/loci/{loc.id}/view?w=600&h=600&pid={display}&reads=hover").json()["background"]
                check(bg["floor"] == rec_pid and bg["cells"] == 1,
                      f"display pid {display}: floor {bg['floor']} (recruitment {rec_pid}), {bg['cells']} cells")


def main():
    tests = [test_canonical_kmers_match_their_reverse_complement, test_diagonal_and_repeat_copies,
             test_hash_sampling_keeps_the_picture, test_tandem_array_is_not_read_as_identical,
             test_uniform_repeat_bin_does_not_light_a_row, test_view_draws_from_the_recruitment_pid,
             test_image_shows_cells_from_the_floor_up]
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
