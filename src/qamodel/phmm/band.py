"""Which cells of the DP to compute.

A band is nothing but two arrays, `lo[u]` and `hi[u]`, giving the half-open
range of template positions evaluated at read position `u`.  Keeping it here,
as pure geometry with no probabilities in sight, means it can be tested against
a known alignment path without running the model at all.

Two constructors, per docs/phmm.md stage F2:

`Band.diagonal` follows the straight diagonal from the origin range.  Width is
the tolerance for cumulative indel drift, so it must be at least as large as the
biggest indel you want representable -- *not* the sequencing-error drift, which
is about two bases over a 15 kb read.  A real 500 bp deletion needs 500 columns
of room no matter how clean the chemistry is.

`Band.from_cigar` follows the aligner's reported path instead, so the band bends
around a large indel and the width only has to cover disagreement with that
path.  Each row spans the whole step the path takes, because a deletion is
crossed within a single row.  This does not reuse the CIGAR to *score* anything (model.tex Remark 6):
it chooses which cells to visit, the score still sums every path inside them,
and `band_edge` reports when the choice was too tight.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

DEFAULT_WIDTH = 100


def _window(arr: np.ndarray, k: int, take) -> np.ndarray:
    """Running min or max over a +/-k window, edges extended."""
    if k <= 0:
        return arr
    return take(sliding_window_view(np.pad(arr, k, mode="edge"), 2 * k + 1),
                axis=1)

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

# S and H are treated as diagonal continuation rather than as consuming nothing.
# A soft-clipped read start is not evidence that those bases came from nowhere;
# it is the aligner declining to place them, and their template almost always
# sits immediately before the reported start. Holding v still across a 500 bp
# clip would put the true diagonal outside any sane width.
_CONSUMES = {
    "M": (1, 1), "=": (1, 1), "X": (1, 1), "S": (1, 1), "H": (1, 1),
    "I": (1, 0),
    "D": (0, 1), "N": (0, 1),
    "P": (0, 0),
}


@dataclass(frozen=True)
class Band:
    """Half-open template range per read position; `lo[u] <= v < hi[u]`."""

    lo: np.ndarray
    hi: np.ndarray

    def __post_init__(self) -> None:
        if self.lo.shape != self.hi.shape:
            raise ValueError("lo and hi must have the same length")
        if np.any(self.hi <= self.lo):
            bad = int(np.argmax(self.hi <= self.lo))
            raise ValueError(f"empty band at read position {bad}")

    @property
    def rows(self) -> int:
        return self.lo.shape[0]

    @property
    def max_width(self) -> int:
        return int((self.hi - self.lo).max())

    def covers(self, u: int, v: int) -> bool:
        return bool(self.lo[u] <= v < self.hi[u])

    @staticmethod
    def _clip(lo, hi, T):
        lo = np.clip(lo, 0, T)
        hi = np.clip(hi, 0, T + 1)
        hi = np.maximum(hi, lo + 1)
        return Band(lo.astype(np.int64), hi.astype(np.int64))

    @staticmethod
    def diagonal(L: int, T: int, origin_lo: int = 0,
                 origin_hi: int | None = None,
                 width: int = DEFAULT_WIDTH) -> "Band":
        """Straight diagonals from every origin in `[origin_lo, origin_hi)`.

        The band spans the whole origin range as well as the drift tolerance, so
        its width is `(origin_hi - origin_lo) + 2*width`.  With an unbounded
        origin range that is the whole template and nothing is saved -- banding
        only pays once the origins are pinned near a reported placement.
        """
        origin_hi = T if origin_hi is None else origin_hi
        u = np.arange(L + 1)
        return Band._clip(origin_lo + u - width, origin_hi + u + width + 1, T)

    @staticmethod
    def from_cigar(L: int, T: int, cigar: str, ref_start: int = 0,
                   width: int = DEFAULT_WIDTH) -> "Band":
        """Follow the path a CIGAR describes, then widen by `width` either side.

        `ref_start` is the alignment's template coordinate, in the same frame as
        the template that will be passed to the DP -- so if the window was cut
        out of a contig, subtract the window's own start first.
        """
        path = cigar_path(L, cigar, ref_start)
        # A deletion is walked entirely inside one row: u does not advance while
        # template is consumed, so the D scan has to cross the whole gap without
        # ever leaving that row's band. Spanning consecutive path points is what
        # makes that possible. Without it a 300 bp deletion leaves rows 499 and
        # 500 with non-overlapping bands and the score collapses to -inf.
        prev = np.concatenate([path[:1], path[:-1]])
        lo = np.minimum(path, prev)
        hi = np.maximum(path, prev)
        # Then spread each row's reach over +/-width rows. Without this the gap
        # is crossable only on the exact row the CIGAR chose, so the equally good
        # paths that place the same event a row earlier or later are cut -- worth
        # log 2 on a two-way tie, and silently, since nothing touches an edge.
        lo = _window(lo, width, np.min)
        hi = _window(hi, width, np.max)
        return Band._clip(lo - width, hi + width + 1, T)


def cigar_path(L: int, cigar: str, ref_start: int = 0) -> np.ndarray:
    """Template position reached after each read base: `v[u]`, length `L+1`.

    Deletions advance `v` without advancing `u`, so they are folded into the
    position recorded for the read base that follows.
    """
    ops = _CIGAR.findall(cigar)
    if not ops or "".join(f"{n}{o}" for n, o in ops) != cigar.strip():
        raise ValueError(f"unparsable CIGAR: {cigar!r}")

    # A leading clip must back the path off the reported start, not run forward
    # from it: with `5S5M` at ref_start 100 the first matched base belongs at
    # template 100, so the five clipped bases occupy 95..99.
    lead = 0
    for count, op in ops:
        if op in ("S", "H"):
            lead += int(count)
        else:
            break

    path = np.empty(L + 1, dtype=np.int64)
    u, v = 0, ref_start - lead
    path[0] = v
    for count, op in ops:
        dq, dr = _CONSUMES[op]
        count = int(count)
        if dq and dr:
            n = min(count, L - u)
            path[u + 1:u + 1 + n] = v + 1 + np.arange(n)
            u += n
            v += count
        elif dq:
            n = min(count, L - u)
            path[u + 1:u + 1 + n] = v
            u += n
        elif dr:
            v += count
            if u <= L:
                path[u] = v          # deletion lands before the next read base
    if u < L:
        # CIGAR shorter than the read: extend the diagonal rather than pinning.
        path[u + 1:] = v + 1 + np.arange(L - u)
    return np.clip(path, 0, None)
