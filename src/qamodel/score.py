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

import hashlib
import math
import os
import sqlite3
from dataclasses import dataclass, field

import numpy as np

from qamodel.exposure import Exposure, Library, exposure_term
from qamodel.phmm.params import PhmmParams, QualityTable, quality_table
from qamodel.phmm.windows import Placement, build_windows, straddles, touches_unit

DEFAULT_C = 200.0
NEG_INF = -math.inf

# A placement whose pid * plen trails the read's best in the same candidate by
# more than this is not scored -- a speed-up only: the boundary rule still sees
# every placement that passes the thresholds, so the set of reads scored does not
# change. pid is the
# aligner's, which weighs every edit alike where the pair HMM weighs them by
# quality and context, so the two disagree by ~10-15 edits' worth. At 0.001 that
# moved unclamped log-ratios by up to 24 nats on a repeat-family locus; at 0.005
# (50 edits on a 10 kb read) it changed no read's ratio by 0.01 while skipping
# 40% of placements (h1tg000067l, 400 reads). An empirical margin, not a bound.
DEFAULT_PRUNE = 0.005
MAX_WORKERS = 32


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
    placements: int = 0                                  # offered, both candidates
    pruned: int = 0                                      # skipped by the pid * plen rule
    boundary: list[str] = field(default_factory=list)    # boundary reads among `reads`
    # boundary reads not scored: no placement touches one candidate's segments
    dropped_boundary: list[str] = field(default_factory=list)

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
            f"boundary reads scored {len(self.boundary)}   dropped (out of one "
            f"candidate's segments) {len(self.dropped_boundary)}",
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
            f"  placements scored {self.placements - self.pruned:,} of {self.placements:,}"
            f" ({self.pruned:,} trailed the read's best pid*plen and were skipped)",
            "",
            "  sensitivity to C:  " + "   ".join(
                f"C={c:g}: {self.at(c):+.1f}" for c in (50, 100, 200, 500, 1000)),
        ]
        return "\n".join(lines)


def _window_key(w) -> bytes:
    h = hashlib.sha1(w.template.encode())
    h.update(np.asarray(w.origin, dtype=np.int64).tobytes())
    h.update(w.band.lo.tobytes())
    h.update(w.band.hi.tobytes())
    return h.digest()


def score_read(read, qual, windows, params: PhmmParams, qt: QualityTable,
               engine, cache: dict | None = None) -> tuple[float, int, float]:
    """log sum over a read's windows in one candidate, window count, worst band edge.

    `cache`, shared across both candidates for one read, skips windows whose
    template, origin and band repeat one already scored -- identical copies of a
    repeat give identical forward sums.
    """
    if not windows:
        return NEG_INF, 0, 0.0
    lps, worst = [], 0.0
    for w in windows:
        key = _window_key(w) if cache is not None else None
        if key is not None and key in cache:
            logp, edge = cache[key]
        else:
            s = engine(read, qual, w.template, params, qt, origin=w.origin, band=w.band)
            logp, edge = s.logp, s.band_edge
            if key is not None:
                cache[key] = (logp, edge)
        lps.append(logp)
        worst = max(worst, edge)
    return _lse(lps), len(windows), worst


def placements(db: sqlite3.Connection, assembly: str, read_name: str,
               prune: float | None = None, pid_min: float | None = None,
               plen_min: float | None = None) -> tuple[list[Placement], list[Placement]]:
    """A read's alignments in one candidate: all that count, and those worth scoring.

    `pid_min` and `plen_min` are the recruitment's thresholds: an alignment the
    recruiter would not follow is not a placement at all -- not scored, and not
    able to drop the read by crossing an edge. A candidate where none passes
    gives the read no window there, so it fits the other one only. With `prune`,
    placements whose pid * plen trails the read's best in this candidate by more
    than `prune` are left out of the second list (see DEFAULT_PRUNE).
    """
    rows = db.execute(
        "SELECT contig, start, end, strand, cigar, soft_left + hard_left, "
        "soft_right + hard_right, pid, plen FROM aln "
        "WHERE assembly = ? AND read_name = ? AND contig IS NOT NULL AND cigar IS NOT NULL",
        (assembly, read_name)).fetchall()
    if pid_min is not None or plen_min is not None:
        rows = [r for r in rows
                if (pid_min is None or (r[7] is not None and r[7] >= pid_min))
                and (plen_min is None or (r[8] is not None and r[8] >= plen_min))]
    every = [Placement(*r[:7]) for r in rows]
    if prune is None or not rows:
        return every, every
    q = [(r[7] or 0.0) * (r[8] or 0.0) for r in rows]
    floor = max(q) - prune
    return every, [p for p, qi in zip(every, q) if qi >= floor]


def default_workers() -> int:
    """QAMODEL_SCORE_WORKERS if set, else the usable cores, capped at MAX_WORKERS."""
    env = os.environ.get("QAMODEL_SCORE_WORKERS")
    if env:
        return max(1, int(env))
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    return max(1, min(MAX_WORKERS, n))


def _default_engine():
    from qamodel.phmm.kernel import HAVE_NUMBA, forward_numba
    from qamodel.phmm.dp import forward_window
    return forward_numba if HAVE_NUMBA else forward_window


def _score_names(names, segs, asm1, asm2, db_path, bams, fastas, params, qt,
                 k, extra, width, engine, prune, pid_min=None, plen_min=None,
                 boundary=frozenset(), progress=None, offset=0, total=None):
    """Score `names` serially; returns a dict of the LocusScore lists and counts:
    scored, dropped, unscorable, boundary, dropped_boundary, offered, pruned.

    Opens its own database, FASTA and read-store handles, so it runs unchanged in
    a worker process.
    """
    import pysam

    from qamodel.phmm.reads import ReadStore

    engine = engine or _default_engine()
    total = len(names) if total is None else total
    db = sqlite3.connect(db_path)
    fa = {a: pysam.FastaFile(fastas[a]) for a in (asm1, asm2)}
    out = {"scored": [], "dropped": [], "unscorable": [], "boundary": [],
           "dropped_boundary": [], "offered": 0, "pruned": 0}
    try:
        with ReadStore(db_path, {a: bams[a] for a in (asm1, asm2)}) as store:
            for i, name in enumerate(names):
                seq, qual = store.get(name)
                pl = {a: placements(db, a, name, prune, pid_min, plen_min) for a in (asm1, asm2)}
                cross = name in boundary
                if cross:
                    # A boundary read crosses a sealed edge, past which both
                    # candidates agree; compared only where it reaches both units.
                    if not all(touches_unit(pl[a][0], segs[a]) for a in (asm1, asm2)):
                        out["dropped_boundary"].append(name)
                        continue
                elif any(straddles(pl[a][0], segs[a]) for a in (asm1, asm2)):
                    out["dropped"].append(name)
                    continue
                ws = {}
                for a in (asm1, asm2):
                    every, kept = pl[a]
                    out["offered"] += len(every)
                    out["pruned"] += len(every) - len(kept)
                    ws[a] = build_windows(kept, segs[a], len(seq), fa[a],
                                          k=k, extra=extra, width=width, cross_edges=cross)
                cache: dict = {}
                lp1, n1, e1 = score_read(seq, qual, ws[asm1].windows, params, qt, engine, cache)
                lp2, n2, e2 = score_read(seq, qual, ws[asm2].windows, params, qt, engine, cache)
                if lp1 == NEG_INF and lp2 == NEG_INF:
                    out["unscorable"].append(name)
                    continue
                out["scored"].append(ReadScore(name, len(seq), lp1, lp2, n1, n2, max(e1, e2)))
                if cross:
                    out["boundary"].append(name)
                if progress:
                    progress(offset + i + 1, total)
    finally:
        for f in fa.values():
            f.close()
        db.close()
    return out


def score_locus(read_names, segments, asm1: str, asm2: str, *,
                db_path: str, bams: dict[str, str], fastas: dict[str, str],
                library: Library, C: float = DEFAULT_C,
                params: PhmmParams | None = None, qt: QualityTable | None = None,
                k: int = 20, extra: int = 50, width: int = 100,
                engine=None, progress=None, prune: float | None = DEFAULT_PRUNE,
                workers: int | None = None, boundary_reads=(),
                pid_threshold: float | None = None,
                plen_threshold: float | None = None) -> LocusScore:
    """Score a recruited locus: read names and segments, as the recruiter returns them.

    `segments` may be any objects with `.assembly`, `.contig`, `.start`, `.end`,
    so this takes qatool Segments without importing qatool. The same read set is
    scored against both candidates; a read crossing a segment edge in either is
    dropped from both, which is what keeps the containment rule symmetric.

    `boundary_reads` are the recruiter's reads crossing a sealed edge: scored too,
    with windows past the edge (the flank is shared, so it cancels from each
    ratio and the exposure term is unchanged), unless one candidate's segments are
    out of their reach. `pid_threshold` and `plen_threshold` are the recruitment's:
    only alignments passing them are placements (see `placements`); None keeps
    every alignment, as before they existed.

    Reads are independent, so with `workers` > 1 they are split into chunks and
    scored in worker processes; results are reassembled in read-name order and
    are identical to a serial run. `workers=None` uses `default_workers()`.
    `prune=None` scores every placement.
    """
    params = params or PhmmParams()
    qt = qt or quality_table()
    workers = default_workers() if workers is None else max(1, workers)

    segs = {asm1: [s for s in segments if s.assembly == asm1],
            asm2: [s for s in segments if s.assembly == asm2]}
    exposure = exposure_term(segs[asm1], segs[asm2], library)
    boundary = frozenset(boundary_reads) - set(read_names)
    names = sorted(set(read_names) | boundary)
    common = (segs, asm1, asm2, db_path, bams, fastas, params, qt, k, extra, width,
              engine, prune, pid_threshold, plen_threshold, boundary)

    # Small chunks balance the load: a read in a repeat family has ~40 windows,
    # a unique read two.
    chunk = max(4, math.ceil(len(names) / (workers * 8)))
    if workers == 1 or len(names) <= 2 * chunk:
        parts = [_score_names(names, *common, progress=progress)]
    else:
        parts = _score_parallel(names, common, workers, chunk, progress)

    tot = {k_: [] for k_ in ("scored", "dropped", "unscorable", "boundary", "dropped_boundary")}
    offered = skipped = 0
    for part in parts:
        for k_ in tot:
            tot[k_] += part[k_]
        offered += part["offered"]
        skipped += part["pruned"]
    return LocusScore(asm1, asm2, tot["scored"], exposure, C, tot["dropped"], tot["unscorable"],
                      placements=offered, pruned=skipped, boundary=tot["boundary"],
                      dropped_boundary=tot["dropped_boundary"])


def _score_parallel(names, common, workers, chunk, progress):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed

    chunks = [names[i:i + chunk] for i in range(0, len(names), chunk)]
    parts = [None] * len(chunks)
    done = 0
    # fork: workers inherit the imported modules and any compiled kernel, and
    # every handle they use is opened inside _score_names.
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=min(workers, len(chunks)), mp_context=ctx) as pool:
        futs = {pool.submit(_score_names, c, *common): j for j, c in enumerate(chunks)}
        for f in as_completed(futs):
            j = futs[f]
            parts[j] = f.result()
            done += len(chunks[j])
            if progress:
                progress(done, len(names))
    return parts
