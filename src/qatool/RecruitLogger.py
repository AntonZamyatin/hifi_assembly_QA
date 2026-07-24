# recruit_logger.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence
import json
import time


@dataclass(frozen=True)
class SupportAln:
    """
    Minimal alignment params you asked to store:
      (read_name, aln_type, pid, plen)
    q_aln_len optional but useful; keep it if you want.
    """
    read_name: str
    aln_type: str
    pid: float
    plen: Optional[float] = None


class RecruitLogger:
    """
    Two TSV files:

    1) segments_added.tsv
       Each line = a delta-segment that was actually added to SegmentStorage,
       + parent segment that expanded into it,
       + JSON list of supporting alignments (read/type/pid/plen).

    2) expanded_alignments.tsv
       Each line = one alignment that passed your filters during expansion from a segment.
    """

    def __init__(self, segments_added_path: str, expanded_alignments_path: str):
        self._seg = open(segments_added_path, "w", buffering=1)
        self._aln = open(expanded_alignments_path, "w", buffering=1)

        self._seg.write(
            "ts\tseg_asm\tseg_contig\tseg_start\tseg_end\t"
            "from_asm\tfrom_contig\tfrom_start\tfrom_end\t"
            "support_json\n"
        )
        self._aln.write(
            "ts\tfrom_asm\tfrom_contig\tfrom_start\tfrom_end\t"
            "aln_contig\taln_start\taln_end\t"
            "read_name\taln_type\tpid\tplen\n"
        )

    def close(self) -> None:
        try:
            self._seg.close()
        except Exception:
            pass
        try:
            self._aln.close()
        except Exception:
            pass

    # -------- File 2: alignment expansion log --------

    def log_expanded_alignment(
        self,
        *,
        from_seg,               # Segment
        aln_contig: str,
        aln_start: int,
        aln_end: int,
        read_name: str,
        aln_type: str,
        pid: float,
        plen: Optional[float],
    ) -> None:
        ts = time.time()
        plen_str = "" if plen is None else f"{float(plen)}"
        self._aln.write(
            f"{ts}\t{from_seg.assembly}\t{from_seg.contig}\t{from_seg.start}\t{from_seg.end}\t"
            f"{aln_contig}\t{int(aln_start)}\t{int(aln_end)}\t"
            f"{read_name}\t{aln_type}\t{float(pid)}\t{plen_str}\n"
        )

    # -------- File 1: added segments log --------

    def log_added_segment(
        self,
        *,
        added_seg,              # Segment (delta that was actually added)
        from_seg,               # Segment that triggered expansion
        support: Sequence[SupportAln],
    ) -> None:
        ts = time.time()
        payload = [
            {
                "name": s.read_name,
                "type": s.aln_type,
                "pid": round(float(s.pid), 6),
                "plen": (None if s.plen is None else round(float(s.plen), 6)),
            }
            for s in support
        ]
        support_json = json.dumps(payload, separators=(",", ":"))
        self._seg.write(
            f"{ts}\t{added_seg.assembly}\t{added_seg.contig}\t{added_seg.start}\t{added_seg.end}\t"
            f"{from_seg.assembly}\t{from_seg.contig}\t{from_seg.start}\t{from_seg.end}\t"
            f"{support_json}\n"
        )