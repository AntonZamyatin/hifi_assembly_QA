"""Sequence and qualities for real reads (stage H of docs/phmm.md).

Secondary alignment records carry no SEQ or QUAL -- the SAM spec allows `*` --
and at a repeat locus the multiplicity reads appear mostly as secondaries. So a
read's sequence always comes from its *primary* record, located through the
`aln` table rather than by scanning the BAM.

SAM stores SEQ in reference orientation. A reverse-strand primary therefore holds
the reverse complement of what was sequenced, with its qualities reversed. The
pair HMM wants the read as sequenced -- the strand belongs to the template, not
the read -- so both are undone here, once, and nowhere else.
"""

from __future__ import annotations

import sqlite3

import numpy as np

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_PRIMARY = ("primary", "primary_wsup")
# Every base's Phred quality for a read stored without any (QUAL '*'), as reads
# aligned from FASTA are -- e.g. homopolymer-compressed ones. The calibration
# anchor of docs/phmm.md stage A: the HiFi data there is mean Q39.
FALLBACK_QUAL = 39


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class ReadStore:
    """`read_name -> (sequence, qualities)`, as sequenced, cached.

    `bams` maps assembly name to BAM path. Any assembly will do, since every read
    is the same read whichever assembly it was aligned to; they are tried in the
    order given, and the first primary record carrying a sequence wins.

    A primary with a sequence but no qualities gets `fallback_qual` on every base
    (None refuses it instead), if no assembly has the read with qualities; such
    reads are collected in `fixed_quality`.
    """

    def __init__(self, db_path: str, bams: dict[str, str],
                 fallback_qual: int | None = FALLBACK_QUAL):
        import pysam

        self._db = sqlite3.connect(db_path)
        self._bams = {a: pysam.AlignmentFile(p, "rb") for a, p in bams.items()}
        self._cache: dict[str, tuple[str, np.ndarray]] = {}
        self.fallback_qual = fallback_qual
        self.fixed_quality: set[str] = set()

    def get(self, read_name: str) -> tuple[str, np.ndarray]:
        hit = self._cache.get(read_name)
        if hit is not None:
            return hit
        no_qual = None
        for asm, bam in self._bams.items():
            row = self._db.execute(
                "SELECT contig, start, strand, q_len FROM aln "
                "WHERE assembly = ? AND read_name = ? AND aln_type IN (?, ?) LIMIT 1",
                (asm, read_name, *_PRIMARY)).fetchone()
            if row is None:
                continue
            contig, start, strand, q_len = row
            for rec in bam.fetch(contig, start, start + 1):
                if (rec.query_name != read_name or rec.is_secondary
                        or rec.is_supplementary or rec.reference_start != start):
                    continue
                seq, qual = rec.query_sequence, rec.query_qualities
                if seq is not None and qual is None and no_qual is None:
                    no_qual = (asm, seq if not rec.is_reverse else revcomp(seq), q_len)
                if seq is None or qual is None:
                    continue
                qual = np.asarray(qual, dtype=np.int64)
                if rec.is_reverse:
                    seq, qual = revcomp(seq), qual[::-1].copy()
                if len(seq) != q_len:
                    # A hard-clipped primary would silently score a truncated read.
                    raise ValueError(
                        f"{read_name}: primary on {asm} stores {len(seq)} bp but the "
                        f"read is {q_len} bp; it is hard-clipped")
                self._cache[read_name] = (seq, qual)
                return seq, qual
        if no_qual is not None:
            # Reads aligned from FASTA (e.g. homopolymer-compressed) carry no QUAL;
            # the pair HMM's emissions need per-base qualities.
            asm, seq, q_len = no_qual
            if self.fallback_qual is None:
                raise ValueError(
                    f"{read_name}: its primary record on {asm} has no base qualities "
                    f"(QUAL '*'), which scoring needs; were the reads aligned from FASTA?")
            if len(seq) != q_len:
                raise ValueError(
                    f"{read_name}: primary on {asm} stores {len(seq)} bp but the "
                    f"read is {q_len} bp; it is hard-clipped")
            qual = np.full(len(seq), self.fallback_qual, dtype=np.int64)
            self.fixed_quality.add(read_name)
            self._cache[read_name] = (seq, qual)
            return seq, qual
        raise KeyError(f"no primary record with a stored sequence for {read_name}")

    def close(self) -> None:
        for bam in self._bams.values():
            bam.close()
        self._db.close()

    def __enter__(self) -> "ReadStore":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
