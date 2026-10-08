"""k-mer similarity between a locus's segments, in bins: the dotplot's background.

For every pair of bins -- asm1 bin `a` (y axis), asm2 bin `b` (x axis) -- the
share of `a`'s distinct canonical k-mers that also occur in `b` estimates their
identity: under substitutions at rate 1 - p a k-mer survives with probability
p^k, so p_hat = share^(1/k). The share is taken both ways -- `a`'s k-mers found
around `b`, and `b`'s found around `a` -- and the smaller one counts. One way
only, a bin of few distinct k-mers is contained in every bin that holds them: a
uniform stretch of a centromeric satellite array lit a whole row against the
other candidate's array (ptg000114l_translocation_001, 928 of 930 columns at
99.7%) while no column was contained in it. Only cells whose estimate reaches a floor are kept,
and the page draws those at or above the locus's recruitment pid -- not the
display pid -- so what shows is sequence close enough for the recruitment to
follow reads across: the diagonal, and every repeat copy that near.

Bins are a regular grid over each axis (plot position // binsize, across segment
boundaries), so the background is one image: it costs the same to draw however
many cells light up -- 878,000 on a 29 Mb locus, which as rectangles would be as
slow as the read lines were. A bin straddling a segment boundary mixes two
segments; with bins far shorter than segments that is a handful of bins.

Three details make the estimate usable:

  The two candidates' grids do not line up: on the diagonal an asm1 bin's k-mers
  fall partly in one asm2 bin and partly in the next. A cell is scored against a
  two-bin window -- its asm2 bin and whichever neighbour gives more -- counting
  each k-mer once, so a split bin still reads as identical; a cell holding less
  than `MIN_SHARE` of its window's matches is not drawn, so the diagonal is not
  doubled. The reverse share uses asm1 windows the same way. Once, not twice: in a repeat the neighbour holds the same k-mers, and
  adding the two bins' counts read a 99.7%-similar satellite array as 100%
  everywhere (h1tg000190l_translocation_001).

  No k-mer is left out for being frequent: in a satellite array every bin is
  near-identical to every other, and that is what the picture should say. The
  cost -- one count per (asm1 bin, asm2 bin) pair holding a shared k-mer, 948
  million on a 4.5 Mb satellite locus -- is bounded instead by hash sampling
  (FracMinHash): past `MAX_PAIRS`, only k-mers whose scrambled code is 0 modulo
  a scale are counted, numerator and denominator alike, which leaves the share
  unbiased; with ~1,000 sampled k-mers in a bin the identity estimate keeps
  about +/-0.0004.

  A bin with fewer than `MIN_KMERS` of a full bin's distinct k-mers is not scored,
  on either axis (mostly N, or the short bin at an axis end).

Cached under `<output>/similarity/`, keyed by the segments, the FASTA files and
the settings.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import pickle

import numpy as np

SIMILARITY_VERSION = 4      # 3: k-mers counted once against a two-bin window; 4: both ways
K = 21
MAX_PAIRS = 200_000_000     # (asm1 bin, asm2 bin) counts per shared k-mer, before sampling
MIN_KMERS = 0.2             # a bin needs this share of a full bin's k-mers to be scored
MIN_SHARE = 0.25            # a cell's share of its (bin, neighbour) pair's matches
FLOOR = 0.95                # cells kept in the cache; the page draws from its own floor up
TARGET_BINS = 1000          # bins on the longer axis, as RegionData's auto binsize
MAX_AXIS_BINS = 4000
_SCRAMBLE = np.uint64(0x9E3779B97F4A7C15)   # odd: a bijection on 64-bit codes

_ENC = np.full(256, 4, dtype=np.uint8)
for _b, _v in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    _ENC[_b] = _v


def auto_binsize(segments, asm1: str, asm2: str) -> int:
    """RegionData's rule: about TARGET_BINS bins on the longer candidate."""
    tot = [sum(s.end - s.start for s in segments if s.assembly == a) for a in (asm1, asm2)]
    return max(1, math.ceil(max(max(tot), 1) / TARGET_BINS))


def axis_map(segments, assembly: str) -> list:
    """(contig, start, end, plot offset) for `assembly`, in the order the plot lays them out."""
    from qatool.plotters.axis import _build_axis_map

    segs = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    return _build_axis_map(segs, assembly)[0]


def axis_length(segs_info) -> int:
    return sum(end - start for _c, start, end, _o in segs_info)


def canonical_kmers(seq: str, k: int = K):
    """(positions, canonical 2-bit codes) of every k-mer free of ambiguous bases."""
    enc = _ENC[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(enc) - k + 1
    if n <= 0:
        return np.empty(0, np.int64), np.empty(0, np.uint64)
    bad = np.concatenate([[0], np.cumsum(enc == 4)])
    ok = (bad[k:] - bad[:-k]) == 0
    e = np.where(enc == 4, 0, enc).astype(np.uint64)
    fwd = np.zeros(n, np.uint64)
    rc = np.zeros(n, np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | e[i:i + n]
        rc |= (np.uint64(3) - e[i:i + n]) << np.uint64(2 * i)
    pos = np.nonzero(ok)[0]
    return pos, np.minimum(fwd, rc)[pos]


def _bin_kmers(fasta, segs_info, binsize: int, k: int):
    """Distinct (k-mer, bin) pairs over a candidate's segments, on the axis's grid,
    sorted by k-mer."""
    keys, bins = [], []
    for contig, start, end, off in segs_info:
        pos, codes = canonical_kmers(fasta.fetch(contig, start, end), k)
        keys.append(codes)
        bins.append(((off + pos) // binsize).astype(np.int64))   # the bin it starts in
    codes = np.concatenate(keys) if keys else np.empty(0, np.uint64)
    b = np.concatenate(bins) if bins else np.empty(0, np.int64)
    order = np.lexsort((b, codes))
    codes, b = codes[order], b[order]
    keep = np.ones(len(codes), bool)
    keep[1:] = (codes[1:] != codes[:-1]) | (b[1:] != b[:-1])
    return codes[keep], b[keep]


def _windows(codes, bins, n: int):
    """Distinct (k-mer, window) pairs, window w = bins w and w + 1: a k-mer in bin j
    belongs to windows j - 1 and j, once each however many of the two bins hold it."""
    w = np.concatenate([bins - 1, bins])
    c = np.concatenate([codes, codes])
    ok = (w >= 0) & (w < n)
    c, w = c[ok], w[ok]
    order = np.lexsort((w, c))
    c, w = c[order], w[order]
    keep = np.ones(len(c), bool)
    keep[1:] = (c[1:] != c[:-1]) | (w[1:] != w[:-1])
    return c[keep], w[keep]


def _match_counts(c1, b1, c2, b2, n1: int, n2: int):
    """M[a, b]: distinct k-mers of asm1 bin a that occur in asm2 bin b."""
    lo = np.searchsorted(c2, c1, side="left")
    hi = np.searchsorted(c2, c1, side="right")
    cnt = hi - lo
    hit = cnt > 0
    rows, lo, cnt = b1[hit], lo[hit], cnt[hit]
    M = np.zeros(n1 * n2, np.int64)
    # in chunks of about `step` expanded pairs, so memory stays bounded
    step = 4_000_000
    cs = np.cumsum(cnt)
    total = int(cs[-1]) if len(cs) else 0
    cuts = np.searchsorted(cs, np.arange(step, total, step), side="left") + 1
    bounds = [0] + sorted({int(c) for c in cuts if 0 < c < len(rows)}) + [len(rows)]
    for i, j in zip(bounds[:-1], bounds[1:]):
        r, l, c = rows[i:j], lo[i:j], cnt[i:j]
        first = np.concatenate([[0], np.cumsum(c)[:-1]])
        idx = np.repeat(l - first, c) + np.arange(c.sum())   # each row's run in c2
        M += np.bincount(np.repeat(r, c) * n2 + b2[idx], minlength=n1 * n2)
    return M.reshape(n1, n2)


def compute(segments, asm1: str, asm2: str, fasta1: str, fasta2: str, *,
            binsize: int | None = None, k: int = K, max_pairs: int = MAX_PAIRS,
            floor: float = FLOOR) -> dict:
    """Cells (asm1 bin a, asm2 bin b) whose estimated identity is at least `floor`.

    Returns the grid (`binsize`, `n1` rows, `n2` columns) and the kept cells as
    index and identity arrays, with their plot edges for convenience.
    """
    import pysam

    segs1, segs2 = axis_map(segments, asm1), axis_map(segments, asm2)
    if not segs1 or not segs2:
        raise ValueError("the locus has no segment in one of the candidates")
    binsize = binsize or auto_binsize(segments, asm1, asm2)
    len1, len2 = axis_length(segs1), axis_length(segs2)
    n1, n2 = math.ceil(len1 / binsize), math.ceil(len2 / binsize)
    if max(n1, n2) > MAX_AXIS_BINS:
        raise ValueError(f"{max(n1, n2)} bins on one axis; raise binsize")
    with pysam.FastaFile(fasta1) as f1, pysam.FastaFile(fasta2) as f2:
        c1, b1 = _bin_kmers(f1, segs1, binsize, k)
        c2, b2 = _bin_kmers(f2, segs2, binsize, k)

    # the cost is one count per (asm1 bin, asm2 bin) pair sharing a k-mer; past
    # max_pairs, count only a hash sample of the k-mers (FracMinHash), both sides alike
    pairs = int((np.searchsorted(c2, c1, side="right") - np.searchsorted(c2, c1, side="left")).sum())
    scale = max(1, math.ceil(pairs / max_pairs))
    if scale > 1:
        keep1 = ((c1 * _SCRAMBLE) % np.uint64(scale)) == 0
        keep2 = ((c2 * _SCRAMBLE) % np.uint64(scale)) == 0
        c1, b1, c2, b2 = c1[keep1], b1[keep1], c2[keep2], b2[keep2]
    kmers1 = np.bincount(b1, minlength=n1)              # distinct (sampled) k-mers per bin
    kmers2 = np.bincount(b2, minlength=n2)
    M = _match_counts(c1, b1, c2, b2, n1, n2)           # k-mers two single bins share
    cw, w = _windows(c2, b2, n2)
    W = _match_counts(c1, b1, cw, w, n1, n2)            # asm1 bin in asm2 two-bin windows
    cw, w = _windows(c1, b1, n1)
    R = _match_counts(c2, b2, cw, w, n2, n1).T          # asm2 bin in asm1 two-bin windows

    # a cell's window: the better of the two that hold the other axis's bin
    # (window w = bins w and w + 1, so w = b - 1 or b)
    prev = np.zeros_like(W)
    prev[:, 1:] = W[:, :-1]
    fwd = np.maximum(W, prev)
    prev = np.zeros_like(R)
    prev[1:, :] = R[:-1, :]
    rev = np.maximum(R, prev)
    full = (binsize - k + 1) / scale
    with np.errstate(divide="ignore", invalid="ignore"):
        share = np.minimum(np.where(kmers1[:, None] > 0, fwd / kmers1[:, None], 0.0),
                           np.where(kmers2[None, :] > 0, rev / kmers2[None, :], 0.0))
        ident = np.where(share > 0, np.power(np.clip(share, 0, 1), 1.0 / k), 0.0)
        own = np.minimum(np.where(fwd > 0, M / fwd, 0.0), np.where(rev > 0, M / rev, 0.0))
    scored = (kmers1 >= MIN_KMERS * full)[:, None] & (kmers2 >= MIN_KMERS * full)[None, :]
    keep = scored & (ident >= floor) & (own >= MIN_SHARE)
    a, b = np.nonzero(keep)
    return {
        "binsize": binsize, "n1": n1, "n2": n2, "len1": len1, "len2": len2,
        "a": a.astype(np.int32), "b": b.astype(np.int32),
        "identity": ident[a, b].astype(np.float32),
        "x0": b * binsize, "x1": np.minimum((b + 1) * binsize, len2),
        "y0": a * binsize, "y1": np.minimum((a + 1) * binsize, len1),
        "settings": {"k": k, "binsize": binsize, "floor": floor, "min_share": MIN_SHARE,
                     "sampling_scale": scale, "pairs": pairs},
    }


def image_rgba(res: dict, floor: float, colors=("#c6dbef", "#08306b")) -> np.ndarray:
    """The kept cells at or above `floor` as an (n1, n2) RGBA image (uint32, row 0
    at the bottom, as Bokeh draws it): light at `floor`, dark at 100%."""
    lo, hi = (np.array([int(c[i:i + 2], 16) for i in (1, 3, 5)], float) for c in colors)
    img = np.zeros((res["n1"], res["n2"]), dtype=np.uint32)
    sel = res["identity"] >= floor
    if sel.any():
        t = (res["identity"][sel] - floor) / max(1e-9, 1.0 - floor)
        t = np.clip(t, 0.0, 1.0)[:, None]
        rgb = np.round(lo + (hi - lo) * t).astype(np.uint32)
        # little-endian RGBA packed into uint32: R | G << 8 | B << 16 | A << 24
        img[res["a"][sel], res["b"][sel]] = (rgb[:, 0] | (rgb[:, 1] << 8) | (rgb[:, 2] << 16)
                                             | (np.uint32(255) << 24))
    return img


def cache_key(segments, asm1: str, asm2: str, fasta1: str, fasta2: str, **settings) -> str:
    payload = json.dumps({
        "version": SIMILARITY_VERSION, "asm": [asm1, asm2], "fasta": [fasta1, fasta2],
        "segments": sorted([s.assembly, s.contig, int(s.start), int(s.end)] for s in segments),
        "k": K, "max_pairs": MAX_PAIRS, "min_kmers": MIN_KMERS, "min_share": MIN_SHARE,
        "floor": FLOOR, **settings,
    }, sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def cache_path(cache_dir: str, *args, **kwargs) -> str:
    return os.path.join(cache_dir, f"similarity_{cache_key(*args, **kwargs)}.pkl")


def load_or_none(path: str) -> dict | None:
    try:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    except (OSError, EOFError, pickle.UnpicklingError):
        return None


def save(path: str, result: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "wb") as fh:
        pickle.dump(result, fh)
    os.replace(tmp, path)
