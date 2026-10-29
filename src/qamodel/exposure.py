"""The count half of the score: how many reads the claimed sequence should draw.

model.tex §7. Each candidate claims some amount of sequence, and claimed sequence
has to pay for itself in reads. The term is `-2 * lambda * L_eff` per candidate,
so the difference between two candidates is

    exposure = -2 * lambda * (L_eff_1 - L_eff_2)

in nats, positive in favour of the first. Without it the emission term alone
prefers a false duplication by exactly log 2 per read (Proposition 6), because
every read gains a second identical origin in the fabricated copy.

lambda: read starts per base per strand, estimated genome-wide and never from
the locus under test -- a locally fitted rate profiles away the very contrast
that identifies copy number (Proposition 3). Every read starts somewhere on the
diploid genome, on one of two strands, so N reads over a diploid length L_dip
give `lambda = N / (2 * L_dip)`. On HG002 that is 1.05e-3, 32x per haplotype.

L_eff: the containment rule scores a read only if it lies wholly inside one
segment, so a segment of length L admits `L - l + 1` starts for a read of span
l, and none if l > L. Averaged over the read-length distribution f and summed
over segments (Proposition 8, exact form):

    L_eff = sum_i sum_l f(l) * max(0, L_i - l + 1)

Read length stands in for reference span, which differs by the net indel count;
that residual belongs to eta_exposure and is small at HiFi error rates.
"""

from __future__ import annotations

import json
import os
import sqlite3
from dataclasses import dataclass

import numpy as np

_PRIMARY = ("primary", "primary_wsup")
_MAX_BIN = 5000                 # 500 kb at 100 bp bins; longer reads are pooled


@dataclass(frozen=True)
class Library:
    """What the sequencing run looks like, genome-wide."""

    n_reads: int                # mapped primary reads
    genome_length: int          # diploid length lambda is measured against
    bin_bp: int
    hist: np.ndarray            # read counts per length bin

    @property
    def lam(self) -> float:
        return self.n_reads / (2.0 * self.genome_length)

    @property
    def centres(self) -> np.ndarray:
        return np.arange(self.hist.size) * self.bin_bp + self.bin_bp / 2

    @property
    def freq(self) -> np.ndarray:
        return self.hist / self.hist.sum()

    @property
    def mean_length(self) -> float:
        return float((self.freq * self.centres).sum())

    @property
    def coverage(self) -> float:
        """Per-haplotype base coverage, for sanity: c = 2 * lambda * mean length."""
        return 2.0 * self.lam * self.mean_length


def fasta_length(path: str) -> int:
    """Total sequence length from a FASTA's .fai index."""
    with open(path + ".fai") as fh:
        return sum(int(line.split("\t")[1]) for line in fh)


def library_stats(db_path: str, assembly: str, genome_length: int, *,
                  bin_bp: int = 100, cache: bool = True) -> Library:
    """Count primary reads and histogram their lengths; cached beside the DB.

    The scan reads every primary row for `assembly` -- about 16 s for HG002's
    12.6 million reads -- so the counts are kept in `<db>.library.json`.
    """
    path = db_path + ".library.json"
    key = f"{assembly}/{bin_bp}"
    store = {}
    if cache and os.path.exists(path):
        with open(path) as fh:
            store = json.load(fh)
        if key in store:
            e = store[key]
            return Library(e["n_reads"], genome_length, bin_bp,
                           np.asarray(e["hist"], dtype=np.int64))

    db = sqlite3.connect(db_path)
    cur = db.execute("SELECT q_len FROM aln WHERE assembly = ? AND aln_type IN (?, ?)",
                     (assembly, *_PRIMARY))
    n, hist = 0, np.zeros(_MAX_BIN, dtype=np.int64)
    while True:
        chunk = cur.fetchmany(500_000)
        if not chunk:
            break
        q = np.fromiter((r[0] for r in chunk), dtype=np.int64, count=len(chunk))
        n += q.size
        hist += np.bincount(np.minimum(q // bin_bp, _MAX_BIN - 1), minlength=_MAX_BIN)
    db.close()
    if n == 0:
        raise ValueError(f"no primary reads for assembly {assembly!r} in {db_path}")

    last = int(np.flatnonzero(hist).max()) + 1
    hist = hist[:last]
    if cache:
        store[key] = {"n_reads": int(n), "hist": hist.tolist()}
        try:
            with open(path, "w") as fh:
                json.dump(store, fh)
        except OSError:                 # a DB shared read-only: recount next time
            pass
    return Library(int(n), genome_length, bin_bp, hist)


def effective_length(lengths, lib: Library) -> float:
    """Starts that keep a read wholly inside its segment, summed over segments."""
    c, f = lib.centres, lib.freq
    return float(sum((f * np.maximum(0.0, L - c + 1.0)).sum() for L in lengths))


@dataclass(frozen=True)
class Exposure:
    L1: int
    L2: int
    k1: int                     # segment counts: fragmentation does not cancel
    k2: int
    Leff1: float
    Leff2: float
    lam: float

    @property
    def value(self) -> float:
        """-2 lambda (L_eff_1 - L_eff_2); positive favours candidate 1."""
        return -2.0 * self.lam * (self.Leff1 - self.Leff2)


def exposure_term(segments1, segments2, lib: Library) -> Exposure:
    l1 = [s.end - s.start for s in segments1]
    l2 = [s.end - s.start for s in segments2]
    return Exposure(sum(l1), sum(l2), len(l1), len(l2),
                    effective_length(l1, lib), effective_length(l2, lib), lib.lam)
