"""Panel application: the three parameter tiers over the Bokeh report.

The tiers exist because the underlying operations differ in cost by four orders of
magnitude, and pretending otherwise is what makes an interactive tool feel broken:

| tier        | cost      | trigger                              |
|-------------|-----------|--------------------------------------|
| recruitment | minutes   | explicit button, subprocess + log    |
| fetch/class | ~1 s      | `value_throttled` (fires on release) |
| display     | instant   | immediate, reuses the cached RegionData |

Display changes never re-read a BAM; Tier-2 changes rebuild `RegionData` but never
re-recruit; recruitment results are cached on disk by parameter hash.

Run:  panel serve src/qaview/app.py --show
  or  python -m qaview --port 5006
"""

import os
import time

import panel as pn
from bokeh.models import GridPlot

from qaview.glyphs import report_layout
from qaview.recruit import RecruitJob, cache_path
from qaview.region import (
    LAYER_COMMON, LAYER_OTHER, LAYER_SUPPLEMENTARY, LAYER_SUPPORT, RegionData,
)
from qaview.state import ExportJob, ViewState

pn.extension(notifications=True)

LAYER_LABELS = {
    LAYER_COMMON: "common",
    LAYER_OTHER: "unclassified",
    LAYER_SUPPORT: "support",
    LAYER_SUPPLEMENTARY: "supplementary",
}

DEFAULTS = dict(
    base="examples", label="hg002", asm1="hifiasm_hic", asm2="curated",
    contig="h2tg000020l", start=61889136, end=61929137,
    recruit_pid=0.997, pid=0.997, plen=0.999,
    best_pid_delta=0.001, pid_tie_eps=0.0001,
    merge_gap=50_000, ext_window=100, ext_max=300_000, ext_rounds=100, max_rounds=100,
)


def dataset_paths(base: str, label: str, asm1: str, asm2: str) -> dict:
    return dict(
        db=f"{base}/dbs/{label}.aln.sqlite",
        bam1=f"{base}/alignments/{label}.{asm1}.sorted.bam",
        bam2=f"{base}/alignments/{label}.{asm2}.sorted.bam",
        fa1=f"{base}/bgz/{label}.{asm1}.fasta.bgz",
        fa2=f"{base}/bgz/{label}.{asm2}.fasta.bgz",
    )


class QaviewApp:
    """Holds the recruitment result and the derived RegionData between interactions."""

    def __init__(self, cache_dir: str = ".qaview_cache", **overrides):
        cfg = {**DEFAULTS, **overrides}
        self.cache_dir = cache_dir
        self.recruitment = None      # {"segments": [...], "reads": {...}}
        self.rd: RegionData | None = None
        self.job: RecruitJob | None = None
        self._poller = None

        w = pn.widgets
        # ---- Tier 1: recruitment -------------------------------------------------
        self.base = w.TextInput(name="base dir", value=cfg["base"])
        self.label = w.TextInput(name="label", value=cfg["label"])
        self.asm1 = w.TextInput(name="asm1 (Y)", value=cfg["asm1"])
        self.asm2 = w.TextInput(name="asm2 (X)", value=cfg["asm2"])
        self.contig = w.TextInput(name="seed contig", value=cfg["contig"])
        self.start = w.IntInput(name="start", value=cfg["start"], step=1000)
        self.end = w.IntInput(name="end", value=cfg["end"], step=1000)
        self.recruit_pid = w.FloatInput(name="recruit pid", value=cfg["recruit_pid"], step=0.001)
        self.merge_gap = w.IntInput(name="merge_gap", value=cfg["merge_gap"], step=1000)
        self.ext_window = w.IntInput(name="ext_window", value=cfg["ext_window"], step=10)
        self.ext_max = w.IntInput(name="ext_max", value=cfg["ext_max"], step=10_000)
        self.ext_rounds = w.IntInput(name="ext_rounds", value=cfg["ext_rounds"])
        self.max_rounds = w.IntInput(name="max_rounds", value=cfg["max_rounds"])
        self.recruit_btn = w.Button(name="Recruit", button_type="primary")
        self.cancel_btn = w.Button(name="Cancel", button_type="warning", disabled=True)
        self.busy = pn.indicators.LoadingSpinner(value=False, size=20)
        self.log = w.TextAreaInput(name="recruitment log", height=180,
                                   value="", disabled=True)

        # ---- Tier 2: fetch / classification -------------------------------------
        self.pid = w.FloatSlider(name="display pid", start=0.90, end=1.0,
                                 step=0.0005, value=cfg["pid"])
        self.plen = w.FloatSlider(name="plen", start=0.90, end=1.0,
                                  step=0.0005, value=cfg["plen"])
        self.best_pid_delta = w.FloatInput(name="best_pid_delta",
                                           value=cfg["best_pid_delta"], step=0.0001)
        self.pid_tie_eps = w.FloatInput(name="pid_tie_eps",
                                        value=cfg["pid_tie_eps"], step=0.0001)
        self.binsize = w.IntInput(name="binsize (0 = auto)", value=0, step=100)
        self.tie_warning = pn.pane.Markdown("", styles={"color": "#b2182b"})

        # ---- Tier 3: display -----------------------------------------------------
        self.color_by = w.RadioButtonGroup(name="colour", options=["class", "block"],
                                           value="class")
        self.layers = w.CheckBoxGroup(
            name="layers", value=list(LAYER_LABELS),
            options={v: k for k, v in LAYER_LABELS.items()})

        # ---- Phase 4: export via report_plot ------------------------------------
        self.export_dir = w.TextInput(name="output dir", value=cfg.get("export_dir", "."))
        self.export_fmt = w.Select(name="format", options=["png", "pdf", "svg"], value="png")
        self.export_dpi = w.IntInput(name="dpi", value=300, step=50)
        self.export_view_only = w.Checkbox(name="current zoom window only", value=True)
        self.export_binsize = w.IntInput(name="binsize (0 = auto)", value=0, step=100)
        self.export_btn = w.Button(name="Export figure", button_type="success")
        self.export_log = w.TextAreaInput(name="export log", height=110, value="", disabled=True)

        self.plot = pn.pane.Bokeh(sizing_mode="stretch_width")
        self.status = pn.pane.Markdown("Set a seed region and press **Recruit**.")

        self.recruit_btn.on_click(self._on_recruit)
        self.cancel_btn.on_click(self._on_cancel)
        self.export_btn.on_click(self._on_export)
        self.export_job = None
        self._export_poller = None
        for widget in (self.pid, self.plen):
            widget.param.watch(self._on_tier2, "value_throttled")
        for widget in (self.best_pid_delta, self.pid_tie_eps, self.binsize):
            widget.param.watch(self._on_tier2, "value")
        for widget in (self.color_by, self.layers):
            widget.param.watch(self._on_tier3, "value")

        self._check_tie()

    # ---- parameters ------------------------------------------------------------

    def _params(self) -> dict:
        p = dataset_paths(self.base.value, self.label.value, self.asm1.value, self.asm2.value)
        p.update(
            asm1=self.asm1.value, asm2=self.asm2.value,
            contig=self.contig.value, start=int(self.start.value), end=int(self.end.value),
            pid_threshold=float(self.recruit_pid.value), plen_threshold=float(self.plen.value),
            best_pid_delta=float(self.best_pid_delta.value),
            merge_gap=int(self.merge_gap.value), extension_window=int(self.ext_window.value),
            ext_max=int(self.ext_max.value), ext_rounds=int(self.ext_rounds.value),
            max_rounds=int(self.max_rounds.value),
        )
        return p

    def _check_tie(self):
        """pid_tie_eps >= best_pid_delta makes the discriminating classes unreachable."""
        if float(self.pid_tie_eps.value) >= float(self.best_pid_delta.value):
            self.tie_warning.object = (
                "**pid_tie_eps >= best_pid_delta** — every non-common read will fall "
                "through to `other` and no read can be called asm1_best/asm2_best.")
        else:
            self.tie_warning.object = ""

    # ---- Tier 1 ----------------------------------------------------------------

    def _on_recruit(self, _event=None):
        params = self._params()
        missing = [v for k, v in params.items()
                   if k in ("db", "bam1", "bam2", "fa1", "fa2") and not os.path.exists(v)]
        if missing:
            self.status.object = "**Missing inputs:**\n\n" + "\n".join(f"- `{m}`" for m in missing)
            return

        cached = os.path.exists(cache_path(self.cache_dir, params))
        self.status.object = ("Loading cached recruitment…" if cached
                              else "Recruiting… (minutes; results are cached)")
        self.busy.value = True
        self.recruit_btn.disabled = True
        self.cancel_btn.disabled = cached
        self.job = RecruitJob(params, self.cache_dir).start()
        self._start_polling()

    def _on_cancel(self, _event=None):
        if self.job is not None:
            self.job.cancel()

    def _start_polling(self):
        if self._poller is not None:
            self._poller.stop()
        # Poll rather than push: the job runs in a subprocess, so there is no
        # in-process callback to hook, and polling keeps the server thread free.
        self._poller = pn.state.add_periodic_callback(self._poll, period=500)

    def _poll(self):
        job = self.job
        if job is None:
            return
        self.log.value = job.log
        if not job.done:
            return
        if self._poller is not None:
            self._poller.stop()
            self._poller = None
        self.busy.value = False
        self.recruit_btn.disabled = False
        self.cancel_btn.disabled = True
        if job.error:
            self.status.object = f"**Recruitment failed:** {job.error}"
            return
        if job.result is None:
            self.status.object = "Recruitment cancelled."
            return
        self.recruitment = job.result
        self._rebuild_region()

    # ---- Tier 2 ----------------------------------------------------------------

    def _on_tier2(self, _event=None):
        self._check_tie()
        if self.recruitment is None:
            return
        self._rebuild_region()

    def _rebuild_region(self):
        if self.recruitment is None:
            return
        p = self._params()
        self.busy.value = True
        try:
            self.rd = RegionData.build(
                self.recruitment["segments"], self.recruitment["reads"],
                p["db"], p["bam1"], p["bam2"],
                asm1_fasta=p["fa1"], asm2_fasta=p["fa2"],
                asm1=p["asm1"], asm2=p["asm2"],
                pid_threshold=float(self.pid.value), plen_threshold=float(self.plen.value),
                binsize=(int(self.binsize.value) or None),
                best_pid_delta=float(self.best_pid_delta.value),
                pid_tie_eps=float(self.pid_tie_eps.value),
            )
        except Exception as exc:                                  # surfaced, not swallowed
            self.status.object = f"**Failed to build region:** `{exc}`"
            self.busy.value = False
            return
        self.busy.value = False
        self._render()

    # ---- Tier 3 ----------------------------------------------------------------

    def _on_tier3(self, _event=None):
        if self.rd is not None:
            self._render()

    def _render(self):
        rd = self.rd
        if rd is None:
            return
        self.plot.object = report_layout(
            rd, color_by=self.color_by.value, visible_layers=set(self.layers.value))
        n_seg = len(rd.segments)
        dropped = rd.lines.n_dropped
        total = len(rd.lines) + dropped
        pct = (100.0 * dropped / total) if total else 0.0
        self.status.object = (
            f"**{len(rd.read_names)} reads · {len(rd.lines)} lines · {n_seg} segments · "
            f"binsize {rd.binsize:,}**  \n"
            f"{rd.asm2} (X) {rd.plot_w:,} bp · {rd.asm1} (Y) {rd.plot_h:,} bp"
            + (f"  \n⚠ {dropped} blocks dropped ({pct:.1f}%)" if pct > 25 else ""))

    # ---- Phase 4: export -------------------------------------------------------

    def _main_figure(self):
        """The dotplot inside the current layout, or None before anything is drawn."""
        layout = self.plot.object
        if layout is None:
            return None
        grid = next((c for c in layout.children if isinstance(c, GridPlot)), None)
        if grid is None:
            return None
        return {(r, c): ch for ch, r, c in grid.children}.get((1, 1))

    def _current_ranges(self):
        """Viewport in qaview plot coordinates, as ((xlo, xhi), (ylo, yhi)).

        In a served session Bokeh syncs range edits back to the server model, so this
        reflects what the user actually zoomed to.  Returns (None, None) when the whole
        region is in view, which exports the full extent rather than a rounded copy.
        """
        fig = self._main_figure()
        if fig is None or self.rd is None:
            return None, None
        xr, yr = fig.x_range, fig.y_range
        x = None if (xr.start <= 0 and xr.end >= self.rd.plot_w) else (xr.start, xr.end)
        y = None if (yr.start <= 0 and yr.end >= self.rd.plot_h) else (yr.start, yr.end)
        return x, y

    def view_state(self, clip_to_view: bool = True) -> ViewState | None:
        if self.rd is None or self.recruitment is None:
            return None
        p = self._params()
        x, y = self._current_ranges() if clip_to_view else (None, None)
        region = f"{self.contig.value}:{int(self.start.value):,}-{int(self.end.value):,}"
        return ViewState.from_region(
            self.rd, db=p["db"], bam1=p["bam1"], bam2=p["bam2"],
            fa1=p["fa1"], fa2=p["fa2"],
            read_names=self.recruitment["reads"],
            pid_threshold=float(self.pid.value), plen_threshold=float(self.plen.value),
            best_pid_delta=float(self.best_pid_delta.value),
            pid_tie_eps=float(self.pid_tie_eps.value),
            x_range=x, y_range=y,
            binsize=(int(self.export_binsize.value) or None),
            dpi=int(self.export_dpi.value),
            title=f"{self.label.value}  {region}",
        )

    def _on_export(self, _event=None):
        state = self.view_state(clip_to_view=self.export_view_only.value)
        if state is None:
            self.status.object = "Nothing to export yet — recruit a region first."
            return
        out_dir = self.export_dir.value or "."
        os.makedirs(out_dir, exist_ok=True)
        stamp = time.strftime("%Y%m%d-%H%M%S")
        name = f"{self.label.value}_{self.contig.value}_{stamp}.{self.export_fmt.value}"
        out = os.path.join(out_dir, name)

        self.export_btn.disabled = True
        self.busy.value = True
        self.export_log.value = f"rendering {state.describe()} …"
        self.export_job = ExportJob(state, out, os.path.join(self.cache_dir, "export")).start()

        if self._export_poller is not None:
            self._export_poller.stop()
        self._export_poller = pn.state.add_periodic_callback(self._poll_export, period=500)

    def _poll_export(self):
        job = self.export_job
        if job is None:
            return
        self.export_log.value = job.log
        if not job.done:
            return
        if self._export_poller is not None:
            self._export_poller.stop()
            self._export_poller = None
        self.export_btn.disabled = False
        self.busy.value = False
        if job.error:
            self.status.object = f"**Export failed:** {job.error}"
        else:
            self.status.object = f"Exported `{job.output_path}`"

    # ---- layout ----------------------------------------------------------------

    def sidebar(self):
        return pn.Accordion(
            ("1 · Recruitment (minutes)", pn.Column(
                self.base, self.label, pn.Row(self.asm1, self.asm2),
                self.contig, pn.Row(self.start, self.end),
                self.recruit_pid, self.merge_gap,
                pn.Row(self.ext_window, self.ext_max),
                pn.Row(self.ext_rounds, self.max_rounds),
                pn.Row(self.recruit_btn, self.cancel_btn, self.busy),
                self.log)),
            ("2 · Fetch & classify (~1 s)", pn.Column(
                self.pid, self.plen,
                pn.Row(self.best_pid_delta, self.pid_tie_eps),
                self.binsize, self.tie_warning)),
            ("3 · Display (instant)", pn.Column(
                pn.pane.Markdown("**colour by**"), self.color_by,
                pn.pane.Markdown("**layers**"), self.layers)),
            ("4 · Export (report_plot)", pn.Column(
                pn.pane.Markdown(
                    "Publication figures come from matplotlib `report_plot`, including "
                    "its k-mer dotplot — not from this page."),
                self.export_dir, pn.Row(self.export_fmt, self.export_dpi),
                self.export_binsize, self.export_view_only,
                self.export_btn, self.export_log)),
            active=[0], toggle=False, width=330,
        )

    def view(self):
        return pn.Row(self.sidebar(),
                      pn.Column(self.status, self.plot, sizing_mode="stretch_width"),
                      sizing_mode="stretch_width")


def build_app(**overrides):
    return QaviewApp(**overrides).view()


if pn.state.served:
    build_app().servable(title="qaview")
