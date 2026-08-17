"""ViewState — what is on screen, in genomic terms — and export through `report_plot`.

The interactive view and the published figure use different plot coordinates: qaview
is exact-bp, `report_plot` pads every segment to a whole number of bins (§5 of
docs/qaview.md).  `ViewState` therefore records **genomic** intervals, never plot
coordinates, so each renderer builds its own mapping from the same description.

Export deliberately goes through matplotlib rather than Bokeh: Bokeh's PNG export
needs selenium plus a chromedriver, and its SVG backend mangles text, neither of
which is acceptable for a figure going into a paper.  `report_plot` is already
debugged and emits vector or 300 dpi output.
"""

import json
import os
from dataclasses import asdict, dataclass, field

from qatool.SegmentStorage import Segment, SegmentCollection
from qaview.jobs import SubprocessJob


def visible_segments(segs_info, assembly: str, lo=None, hi=None) -> list:
    """Clip an exact-bp axis map to a plot-coordinate window, as genomic Segments.

    `lo`/`hi` of None mean "no bound on that side", which reproduces the full region.
    Coordinates are rounded rather than truncated: Bokeh range ends arrive as floats
    and a bare int() would shave a base off every segment on every export.
    """
    out = []
    for contig, start, end, off in segs_info:
        seg_len = end - start
        a = off if lo is None else max(lo, off)
        b = (off + seg_len) if hi is None else min(hi, off + seg_len)
        if b <= a:
            continue
        gs = start + int(round(a - off))
        ge = start + int(round(b - off))
        gs = max(start, min(gs, end))
        ge = max(start, min(ge, end))
        if ge > gs:
            out.append(Segment(assembly, contig, gs, ge))
    return out


@dataclass
class ViewState:
    """Everything needed to re-render the current view in either renderer."""
    db: str
    bam1: str
    bam2: str
    fa1: str
    fa2: str
    asm1: str
    asm2: str
    segments: list                      # genomic, both assemblies
    read_names: list
    pid_threshold: float
    plen_threshold: float
    best_pid_delta: float
    pid_tie_eps: float
    classification: str = "support"
    binsize: int | None = None          # None -> report_plot picks, respecting its k floor
    k: int = 15
    max_kmer_freq: int = 100
    support_stripes: bool = True
    dpi: int = 300
    title: str | None = None

    @classmethod
    def from_region(cls, rd, *, db, bam1, bam2, fa1, fa2, read_names,
                    pid_threshold, plen_threshold, best_pid_delta, pid_tie_eps,
                    x_range=None, y_range=None, **kw) -> "ViewState":
        """Snapshot a RegionData, optionally clipped to the current viewport.

        x_range is the asm2 (X) window, y_range the asm1 (Y) window, both in qaview
        plot coordinates; pass None for the whole axis.
        """
        xlo, xhi = x_range if x_range else (None, None)
        ylo, yhi = y_range if y_range else (None, None)
        segments = (visible_segments(rd.asm1_segs_info, rd.asm1, ylo, yhi)
                    + visible_segments(rd.asm2_segs_info, rd.asm2, xlo, xhi))
        return cls(db=db, bam1=bam1, bam2=bam2, fa1=fa1, fa2=fa2,
                   asm1=rd.asm1, asm2=rd.asm2,
                   segments=segments, read_names=sorted(read_names),
                   pid_threshold=pid_threshold, plen_threshold=plen_threshold,
                   best_pid_delta=best_pid_delta, pid_tie_eps=pid_tie_eps, **kw)

    def describe(self) -> str:
        by_asm = {}
        for s in self.segments:
            by_asm.setdefault(s.assembly, []).append(s)
        parts = []
        for asm in (self.asm1, self.asm2):
            segs = by_asm.get(asm, [])
            span = sum(s.end - s.start for s in segs)
            parts.append(f"{asm}: {len(segs)} seg / {span:,} bp")
        return "  ·  ".join(parts)

    def to_json(self) -> str:
        d = asdict(self)
        d["segments"] = [[s.assembly, s.contig, s.start, s.end] for s in self.segments]
        return json.dumps(d)

    @classmethod
    def from_json(cls, text: str) -> "ViewState":
        d = json.loads(text)
        d["segments"] = [Segment(*s) for s in d["segments"]]
        return cls(**d)

    def render(self, output_path: str) -> str:
        """Draw with `report_plot`.  Blocking — run via ExportJob from a server."""
        import matplotlib
        matplotlib.use("Agg")
        from qatool.plotters import report_plot

        report_plot(
            "", "", SegmentCollection(list(self.segments)), set(self.read_names),
            self.db, self.pid_threshold, self.plen_threshold,
            k=self.k, binsize=self.binsize,
            output_path=output_path, dpi=self.dpi,
            asm1_bam=self.bam1, asm2_bam=self.bam2,
            asm1_fasta=self.fa1, asm2_fasta=self.fa2,
            asm1=self.asm1, asm2=self.asm2,
            pid_tie_eps=self.pid_tie_eps, best_pid_delta=self.best_pid_delta,
            max_kmer_freq=self.max_kmer_freq,
            support_stripes=self.support_stripes,
            classification=self.classification,
            title=self.title, show=False, verbose=True,
        )
        return output_path


class ExportJob(SubprocessJob):
    """Render a ViewState to a file in a child process."""

    def __init__(self, state: ViewState, output_path: str, work_dir: str,
                 python: str | None = None):
        self.state = state
        self.output_path = output_path
        os.makedirs(work_dir, exist_ok=True)
        self._state_file = os.path.join(work_dir, "export_state.json")
        super().__init__(["-m", "qaview.export_worker", self._state_file, output_path], python)

    @property
    def label(self) -> str:
        return "export"

    def start(self) -> "ExportJob":
        with open(self._state_file, "w") as fh:
            fh.write(self.state.to_json())
        return super().start()

    def _on_success(self):
        if not os.path.exists(self.output_path):
            raise FileNotFoundError(self.output_path)
