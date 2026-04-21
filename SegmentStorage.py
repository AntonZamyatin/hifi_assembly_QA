from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple
import pysam

@dataclass(frozen=True, order=True)
class Segment:
    assembly: str   # e.g. "asm1" or "asm2"
    contig: str
    start: int
    end: int
    """
    Dataclass for representation of assembly contig substrings
    """

    def __post_init__(self):
        if self.start < 0 or self.end < 0:
            raise ValueError(f"Negative coords: {self}")
        if self.end < self.start:
            raise ValueError(f"Invalid interval (end < start): {self}")

    def __len__(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        return f"({self.assembly} {self.contig}){self.start}:{self.end}"
        
    @property
    def key(self) -> Tuple[str, str]:
        return self.assembly, self.contig

class SegmentCollection:
    """
    Collection of non-overlapping segments that can be added to and iterated over.
    """
    def __init__ (self, segments: List[Segment]|SegmentCollection = None, merge_gap: int = 0):
        self.merge_gap = merge_gap
        self.count = 0
        self.storage: Dict[Tuple[str, str], List[Tuple[int, int]]] = dict()
        if segments:
            for segment in segments:
                self.add(segment)

    def add(self, segment: Segment):
        segment_list: List[Tuple[int, int]] = self.storage.setdefault((segment.assembly, segment.contig), [])
        new_segment_list: List[Tuple[int, int]] = []
        added: bool = False
        seg_start = segment.start
        seg_end = segment.end
        for i, interval in enumerate(segment_list):
            interval_start, interval_end = interval
            if added:
                new_segment_list.append(interval)
                continue
            # completely before interval -> add new segment and interval
            if seg_end < interval_start - self.merge_gap:
                new_segment_list.append((seg_start, seg_end))
                new_segment_list.append(interval)
                added = True
                continue
            # completely after interval -> add interval
            if seg_start >= interval_end + self.merge_gap:
                new_segment_list.append(interval)
                continue
            # overlap -> merge
            seg_start = min(seg_start, interval_start)
            seg_end = max(seg_end, interval_end)
            self.count -= 1
        if not added:
            new_segment_list.append((seg_start, seg_end))
        self.storage[(segment.assembly, segment.contig)] = new_segment_list
        self.count += 1

    def __iter__(self):
        for (asm, contig), intervals in self.storage.items():
            for start, end in intervals:
                yield Segment(asm, contig, start, end)

    def __len__(self):
        return self.count

    def __add__(self, other: SegmentCollection) -> SegmentCollection:
        for segment in other:
            self.add(segment)
        return self

    def __str__(self) -> str:
        out_str = ""
        if not self.is_empty():
            for seg in self:
                out_str += f"{seg}\n"
        return out_str

    def is_empty(self) -> bool:
        return self.count == 0


def get_concatenated_segments(sc: SegmentCollection, ref_fasta: str,
                                    cand_fasta: str, asm1: str = "asm1",
                                    asm2: str = "asm2") -> Tuple[str, str]:
    """
    Concatenate segments in SegmentCollection order into two strings
    (ref, cand) using the provided assembly FASTA files.
    """
    ref_parts: List[str] = []
    cand_parts: List[str] = []

    ref_fa = pysam.FastaFile(ref_fasta)
    cand_fa = pysam.FastaFile(cand_fasta)
    try:
        for seg in sc:
            if seg.assembly == asm1:
                ref_parts.append(ref_fa.fetch(seg.contig, seg.start, seg.end))
            elif seg.assembly == asm2:
                cand_parts.append(cand_fa.fetch(seg.contig, seg.start, seg.end))
            else:
                raise ValueError(f"Unknown assembly label: {seg.assembly}")
    finally:
        ref_fa.close()
        cand_fa.close()

    return "".join(ref_parts), "".join(cand_parts)


if __name__ == "__main__":
    ss = SegmentCollection(merge_gap=1)
    ss.add(Segment("ref", "contig1", 1000, 2000))
    ss.add(Segment("ref", "contig1", 300, 400))
    ss.add(Segment("ref", "contig1", 100, 501))
    ss.add(Segment("ref", "contig3", 500, 600))
    ss.add(Segment("ref", "contig1", 600, 700))
    ss.add(Segment("ref", "contig1", 700, 800))

    for s in ss:
        print(s)
