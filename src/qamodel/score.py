"""The locus score: which of two candidate assemblies the reads support.

    Delta S = sum_j clamp( log sum_{p in A1} e_jp - log sum_{p in A2} e_jp, -C, +C )
              - 2 * lambda * (L_eff_1 - L_eff_2)

Positive favours the first candidate. Units are nats.

The emission term is the pair HMM, summed over every window a read has in each
candidate. The exposure term is `qamodel.exposure`. The clamp bounds any single
read's say at C: without it, a 15 kb read that cannot be placed in one candidate
is explained as insertions at ~2.6 nats a base, about -39,000 nats, and outweighs
five hundred reads that fit both -- by an amount set by its length rather than by
how much sequence is missing. Clamping every read, not only the one-sided ones,
keeps the score continuous: there is no placeability threshold for a read to fall
across. A read that fits only one candidate simply lands on the clamp.

Delta S is not a likelihood ratio once clamped (model.tex Remark 8). It must not
be read as a Bayes factor or turned into posterior odds; its mapping to error
rates is a calibration question for simulation. `LocusScore.at(C)` recomputes it
for another C without rescoring anything, because C is a choice to be swept, not
tuned per locus.
"""

from __future__ import annotations

import math
import sqlite3
from dataclasses import dataclass, field

import numpy as np

from qamodel.exposure import Exposure, Library, exposure_term
from qamodel.phmm.params import PhmmParams, QualityTable, quality_table
from qamodel.phmm.windows import Placement, build_windows

DEFAULT_C = 200.0
NEG_INF = -math.inf


def _lse(xs) -> float:
    xs = [x for x in xs if math.isfinite(x)]
    if not xs:
        return NEG_INF
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def _clamp(x: float, C: float) -> float:
    return max(-C, min(C, x))


@dataclass(frozen=True)
class ReadScore:
    name: str
    length: int
    logp1: float                # log sum over A1's windows; -inf if none
    logp2: float
    windows1: int
    windows2: int
    band_edge: float            # worst over every window scored for this read

    @property
    def ratio(self) -> float:
        """Unclamped log-ratio; +/-inf when the read fits only one candidate."""
        return self.logp1 - self.logp2

    def clamped(self, C: float) -> float:
        return _clamp(self.ratio, C)


@dataclass
class LocusScore:
    asm1: str
    asm2: str
    reads: list[ReadScore]
    exposure: Exposure
    C: float
    dropped: list[str] = field(default_factory=list)     # crossed a segment boundary
    unscorable: list[str] = field(default_factory=list)  # no window in either candidate

    def emission_at(self, C: float) -> float:
        return float(sum(r.clamped(C) for r in self.reads))

    @property
    def emission(self) -> float:
        return self.emission_at(self.C)

    @property
    def delta_s(self) -> float:
        return self.emission + self.exposure.value

    def at(self, C: float) -> float:
        """Delta S for another clamp, from the stored per-read ratios."""
        return self.emission_at(C) + self.exposure.value

    def clamped_counts(self, C: float | None = None) -> tuple[int, int]:
        """Reads sitting on +C and on -C."""
        C = self.C if C is None else C
        return (sum(r.ratio >= C for r in self.reads),
                sum(r.ratio <= -C for r in self.reads))

    @property
    def worst_band_edge(self) -> float:
        return max((r.band_edge for r in self.reads), default=0.0)

    def summary(self) -> str:
        e = x = self.exposure
        up, down = self.clamped_counts()
        favoured = self.asm1 if self.delta_s > 0 else self.asm2
        lines = [
            f"=== pair-HMM locus score: {self.asm1} vs {self.asm2} ===",
            f"reads scored {len(self.reads)}   dropped (crossed a segment edge) "
            f"{len(self.dropped)}   unscorable {len(self.unscorable)}",
            "",
            f"  emission   {self.emission:+12.2f} nats   "
            f"({up} reads at +C, {down} at -C, C = {self.C:g})",
            f"  exposure   {x.value:+12.2f} nats   "
            f"-2 x {e.lam:.3e} x ({e.Leff1:,.0f} - {e.Leff2:,.0f})",
            f"  Delta S    {self.delta_s:+12.2f} nats   -> favours {favoured}",
            "",
            f"  segments   {self.asm1}: {e.k1} totalling {e.L1:,} bp   "
            f"{self.asm2}: {e.k2} totalling {e.L2:,} bp",
            f"  worst band_edge over all windows {self.worst_band_edge:.1e}",
            "",
            "  sensitivity to C:  " + "   ".join(
                f"C={c:g}: {self.at(c):+.1f}" for c in (50, 100, 200, 500, 1000)),
        ]
        return "\n".join(lines)


def score_read(read, qual, windows, params: PhmmParams, qt: QualityTable,
               engine) -> tuple[float, int, float]:
    """log sum over a read's windows in one candidate, window count, worst band edge."""
    if not windows:
        return NEG_INF, 0, 0.0
    lps, worst = [], 0.0
    for w in windows:
        s = engine(read, qual, w.template, params, qt, origin=w.origin, band=w.band)
        lps.append(s.logp)
        worst = max(worst, s.band_edge)
    return _lse(lps), len(windows), worst


def placements(db: sqlite3.Connection, assembly: str, read_name: str) -> list[Placement]:
    rows = db.execute(
        "SELECT contig, start, end, strand, cigar, soft_left + hard_left, "
        "soft_right + hard_right FROM aln "
        "WHERE assembly = ? AND read_name = ? AND contig IS NOT NULL AND cigar IS NOT NULL",
        (assembly, read_name)).fetchall()
    return [Placement(*r) for r in rows]


def score_locus(read_names, segments, asm1: str, asm2: str, *,
                db_path: str, bams: dict[str, str], fastas: dict[str, str],
                library: Library, C: float = DEFAULT_C,
                params: PhmmParams | None = None, qt: QualityTable | None = None,
                k: int = 20, extra: int = 50, width: int = 100,
                engine=None, progress=None) -> LocusScore:
    """Score a recruited locus: read names and segments, as the recruiter returns them.

    `segments` may be any objects with `.assembly`, `.contig`, `.start`, `.end`,
    so this takes qatool Segments without importing qatool. The same read set is
    scored against both candidates; a read crossing a segment edge in either is
    dropped from both, which is what keeps the containment rule symmetric.
    """
    import pysam

    from qamodel.phmm.reads import ReadStore

    if engine is None:
        from qamodel.phmm.kernel import HAVE_NUMBA, forward_numba
        from qamodel.phmm.dp import forward_window
        engine = forward_numba if HAVE_NUMBA else forward_window
    params = params or PhmmParams()
    qt = qt or quality_table()

    segs = {asm1: [s for s in segments if s.assembly == asm1],
            asm2: [s for s in segments if s.assembly == asm2]}
    exposure = exposure_term(segs[asm1], segs[asm2], library)

    db = sqlite3.connect(db_path)
    fa = {a: pysam.FastaFile(fastas[a]) for a in (asm1, asm2)}
    scored, dropped, unscorable = [], [], []
    try:
        with ReadStore(db_path, {a: bams[a] for a in (asm1, asm2)}) as store:
            names = sorted(read_names)
            for i, name in enumerate(names):
                seq, qual = store.get(name)
                ws = {a: build_windows(placements(db, a, name), segs[a], len(seq), fa[a],
                                       k=k, extra=extra, width=width)
                      for a in (asm1, asm2)}
                if ws[asm1].straddles or ws[asm2].straddles:
                    dropped.append(name)
                    continue
                lp1, n1, e1 = score_read(seq, qual, ws[asm1].windows, params, qt, engine)
                lp2, n2, e2 = score_read(seq, qual, ws[asm2].windows, params, qt, engine)
                if lp1 == NEG_INF and lp2 == NEG_INF:
                    unscorable.append(name)
                    continue
                scored.append(ReadScore(name, len(seq), lp1, lp2, n1, n2, max(e1, e2)))
                if progress:
                    progress(i + 1, len(names))
    finally:
        for f in fa.values():
            f.close()
        db.close()

    return LocusScore(asm1, asm2, scored, exposure, C, dropped, unscorable)
