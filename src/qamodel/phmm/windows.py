"""From aligner placements to scoring windows (stage F3 of docs/phmm.md).

A read's origin sum runs over every start position in the recruited unit. The
aligner's placements say where that mass is; this module turns each one into a
window -- a strand-oriented template, a range of starts the window owns, and a
band -- so the pair HMM only visits the handful of places that matter.

Four rules, each of which was found to matter the hard way:

  Containment is judged on the *aligned* span, the same rule the recruiter uses.
  A placement inside a segment is scored; one crossing a segment boundary means
  the read's origin may sit outside the unit, and the whole read is dropped from
  both candidates so the read set stays symmetric. Placements elsewhere in the
  genome are not part of this unit's sum and are ignored.

  Clips are recorded in reference orientation. For a reverse-strand placement the
  read's first base as sequenced is at the reference *right* end, so the read's
  start is `end + clip_right - 1`, and the template is the reverse complement of
  the window with the CIGAR reversed to match.

  Origin ranges must be disjoint, or shared starts are counted twice and every
  repeat is inflated. Placements whose ranges would overlap are merged into one
  window owning the hull. Extra starts inside a hull cost nothing: they explain
  the read poorly and contribute almost no mass.

  The band is the union of each member's CIGAR-anchored band, so it bends around
  large indels the aligner already found, and covers every seeded start.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from qamodel.phmm.band import _CIGAR, Band
from qamodel.phmm.reads import revcomp


def reverse_cigar(cigar: str) -> str:
    """The same alignment read from the other end."""
    return "".join(f"{n}{op}" for n, op in reversed(_CIGAR.findall(cigar)))


@dataclass(frozen=True)
class Placement:
    """One alignment of a read, in reference coordinates, 0-based half-open.

    `clip_left` and `clip_right` are soft plus hard clips in CIGAR order, which
    is reference orientation -- not the read's own left and right.
    """

    contig: str
    start: int
    end: int
    strand: int                 # 1 forward, 0 reverse
    cigar: str
    clip_left: int
    clip_right: int

    @property
    def read_start(self) -> int:
        """Reference position of the read's first base as sequenced."""
        if self.strand == 1:
            return self.start - self.clip_left
        return self.end + self.clip_right - 1


@dataclass(frozen=True)
class Window:
    """One call's worth of work: template, owned starts, and band."""

    contig: str
    strand: int
    ref_lo: int                 # template = assembly[contig][ref_lo:ref_hi],
    ref_hi: int                 # reverse-complemented when strand == 0
    template: str
    origin: tuple[int, int]     # half-open, in template coordinates
    band: Band
    members: int                # placements merged into this window


@dataclass(frozen=True)
class WindowSet:
    windows: list[Window]
    straddles: bool             # a placement crossed a segment boundary
    outside: int                # placements outside the unit, ignored
    lost: int                   # placements whose starts fell off a clipped window


def build_windows(placements, segments, read_len: int, fasta, *,
                  k: int = 20, extra: int = 50, width: int = 100) -> WindowSet:
    """Windows for one read against one candidate.

    `segments` are that candidate's recruited segments -- anything with
    `.contig`, `.start`, `.end`. `fasta` needs only `.fetch(contig, lo, hi)`.
    `k` is how far the true start may sit from where the aligner put it;
    `extra` is template beyond the read's reach, room for deletions; `width` is
    the band's tolerance around the aligner's path.
    """
    groups: dict[tuple, list[Placement]] = {}
    straddles, outside = False, 0
    for p in placements:
        home, touches = None, False
        for s in segments:
            if s.contig != p.contig or p.end <= s.start or p.start >= s.end:
                continue
            touches = True
            if s.start <= p.start and p.end <= s.end:
                home = s
                break
        if home is None:
            if touches:
                straddles = True
            else:
                outside += 1
            continue
        groups.setdefault((home.contig, home.start, home.end, p.strand), []).append(p)

    windows: list[Window] = []
    lost = 0
    for (contig, seg_lo, seg_hi, strand), ps in groups.items():
        ps = sorted(ps, key=lambda p: p.read_start)
        merged = [[ps[0]]]
        for p in ps[1:]:
            # Ranges [s-k, s+k] overlap exactly when starts are within 2k.
            if p.read_start - merged[-1][-1].read_start <= 2 * k:
                merged[-1].append(p)
            else:
                merged.append([p])

        for grp in merged:
            lo = max(seg_lo, min(p.start - p.clip_left for p in grp) - k - extra)
            hi = min(seg_hi, max(p.end + p.clip_right for p in grp) + k + extra)
            tmpl = fasta.fetch(contig, lo, hi).upper()
            if strand == 0:
                tmpl = revcomp(tmpl)
            T = len(tmpl)

            o_lo, o_hi = T, 0
            b_lo = b_hi = None
            for p in grp:
                if strand == 1:
                    o = p.read_start - lo
                    band = Band.from_cigar(read_len, T, p.cigar, p.start - lo, width)
                else:
                    o = hi - 1 - p.read_start
                    band = Band.from_cigar(read_len, T, reverse_cigar(p.cigar),
                                           hi - p.end, width)
                o_lo, o_hi = min(o_lo, o - k), max(o_hi, o + k + 1)
                b_lo = band.lo if b_lo is None else np.minimum(b_lo, band.lo)
                b_hi = band.hi if b_hi is None else np.maximum(b_hi, band.hi)

            o_lo, o_hi = max(o_lo, 0), min(o_hi, T)
            if T == 0 or o_hi <= o_lo:
                lost += len(grp)
                continue
            windows.append(Window(contig, strand, lo, hi, tmpl, (o_lo, o_hi),
                                  Band(b_lo, b_hi), len(grp)))

    return WindowSet(windows, straddles, outside, lost)
