"""Panel application: the three parameter tiers over the Bokeh report.

The tiers exist because the underlying operations differ in cost by four orders of
magnitude, and pretending otherwise is what makes an interactive tool feel broken:

| tier        | cost      | trigger                              |
|-------------|-----------|--------------------------------------|
| recruitment | minutes   | explicit button, subprocess + log    |
| scoring     | ~30 s     | automatic after recruitment, subprocess, cached |
| fetch/class | ~1 s      | `value_throttled` (fires on release) |
| display     | instant   | immediate, reuses the cached RegionData |

The pair-HMM clamp C is a display setting: Delta S at another C comes from the
stored per-read ratios, so the slider never rescores.

Display changes never re-read a BAM; Tier-2 changes rebuild `RegionData` but never
re-recruit; recruitment results are cached on disk by parameter hash.

Run:  panel serve src/qaview/app.py --show
  or  python -m qaview --port 5006
"""

import os
import time

import pandas as pd
import panel as pn
import param
from bokeh.models import GridPlot

from qaview.glyphs import BAR_SIZE, MARGIN_SIZE, RATIO_NEG, RATIO_POS, report_layout
from qaview.recruit import RecruitJob, cache_path
from qaview.score import ScoreJob
from qaview.seeds import load_seeds, scores_dir, seeds_path
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

SIDEBAR_WIDTH = 330
# Pixels around the dotplot that are not the dotplot, excluding the marginals:
# vertically the status/display row, Delta S panel, in-figure readout and the axis
# below; horizontally the sidebar, left contig bar, toolbar and margins.
CHROME_HEIGHT = 280 + BAR_SIZE
CHROME_WIDTH = SIDEBAR_WIDTH + BAR_SIZE + 110
# Marginal size by window height: a full 200 px on tall screens, less on laptops so
# the dotplot keeps most of the height.
SMALL_MARGIN, SMALL_SCREEN = 120, 1100


class Viewport(pn.custom.JSComponent):
    """Reports the browser window size, so the figure can be fitted into it.

    Bokeh figures here have fixed pixel sizes (square pixels on the dotplot matter
    more than filling the box), so the server has to know the box.  Resizes are
    debounced; each one triggers a single redraw.
    """

    width = param.Integer(default=0)
    height = param.Integer(default=0)

    _esm = """
    export function render({ model }) {
      const report = () => {
        model.width = window.innerWidth;
        model.height = window.innerHeight;
      };
      let timer = null;
      window.addEventListener("resize", () => {
        clearTimeout(timer);
        timer = setTimeout(report, 300);
      });
      report();
      return document.createElement("div");
    }
    """


DEFAULTS = dict(
    base="examples", label="hg002", asm1="hifiasm_hic", asm2="curated",
    contig="h2tg000020l", start=61889136, end=61929137,
    recruit_pid=0.997, pid=0.997, plen=0.999,
    best_pid_delta=0.001, pid_tie_eps=0.0001,
    merge_gap=50_000, ext_window=100, ext_max=300_000, ext_rounds=100, max_rounds=100,
    quast_run="hg002_hic",
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
        self.score = None            # qamodel LocusScore for the current recruitment
        self.score_job: ScoreJob | None = None
        self._score_poller = None

        w = pn.widgets
        # ---- Seeds: QUAST misassemblies ------------------------------------------
        self.quast_run = w.TextInput(name="QUAST run (<base>/quast/<run>/misassemblies.tsv)",
                                     value=cfg["quast_run"])
        self.seed_table = w.Tabulator(
            pd.DataFrame(columns=["seed", "bp", "delta_s", "cached"]),
            show_index=False, selectable=1, disabled=True, height=300,
            pagination=None, layout="fit_columns",
            widths={"seed": 128, "bp": 76, "delta_s": 58, "cached": 28},
            formatters={"bp": {"type": "money", "precision": 0, "thousand": ","},
                        "delta_s": {"type": "money", "precision": 1, "thousand": ""}},
            stylesheets=[".tabulator-cell, .tabulator-col-title { font-size: 11px; }"],
            header_filters={"seed": {"type": "input", "func": "like",
                                     "placeholder": "filter"}},
            hidden_columns=["region", "contig", "type", "breakpoint", "start", "end",
                            "inconsistency_bp"],
            titles={"bp": "breakpoint", "delta_s": "ΔS nb", "cached": "✓"})
        self.seed_reload = w.Button(name="reload", width=80)
        self.seed_note = pn.pane.Markdown("", width=SIDEBAR_WIDTH - 140,
                                          styles={"font-size": "12px", "color": "#666"})

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
        # log-ratio by default: before the score lands it draws the support classes.
        self.color_by = w.RadioButtonGroup(name="colour",
                                           options=["class", "block", "log-ratio"],
                                           value="log-ratio", width=260)
        self.clamp = w.DiscreteSlider(name="clamp C (nats per read)",
                                      options=[10, 20, 50, 100, 200, 500, 1000], value=200,
                                      width=260)
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
        self.status = pn.pane.Markdown("Pick a seed or set one, and press **Recruit**.",
                                       sizing_mode="stretch_width")
        self.score_pane = pn.pane.HTML("", sizing_mode="stretch_width")
        self.score_log = w.TextAreaInput(name="scoring log", height=140, value="",
                                         disabled=True)

        self.recruit_btn.on_click(self._on_recruit)
        self.quast_run.param.watch(self._load_seeds, "value")
        self.seed_reload.on_click(self._load_seeds)
        self.base.param.watch(self._load_seeds, "value")
        self.seed_table.param.watch(self._on_seed, "selection")
        self.cancel_btn.on_click(self._on_cancel)
        self.export_btn.on_click(self._on_export)
        self.export_job = None
        self._export_poller = None
        for widget in (self.pid, self.plen):
            widget.param.watch(self._on_tier2, "value_throttled")
        for widget in (self.best_pid_delta, self.pid_tie_eps, self.binsize):
            widget.param.watch(self._on_tier2, "value")
        for widget in (self.color_by, self.layers, self.clamp):
            widget.param.watch(self._on_tier3, "value")
        self.viewport = Viewport(width=0, height=0, margin=0)
        self.viewport.param.watch(self._on_tier3, ["width", "height"])

        self._check_tie()
        self._load_seeds()

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

    # ---- seeds -----------------------------------------------------------------

    def _load_seeds(self, _event=None):
        path = seeds_path(self.base.value, self.quast_run.value)
        if not os.path.exists(path):
            self.seed_table.value = self.seed_table.value.iloc[0:0]
            self.seed_note.object = f"no seed table at `{path}`"
            return
        df = load_seeds(path, scores_dir(self.base.value, self.quast_run.value))
        df["cached"] = self._cached_flags(df)
        self.seed_table.value = df
        n_scored = int(df["delta_s"].notna().sum())
        self.seed_note.object = (
            f"{len(df)} misassemblies on {df['contig'].nunique()} {self.asm1.value} contigs. "
            f"ΔS nb: the notebook's pair-HMM score ({n_scored} so far). "
            "✓: recruitment cached with the current parameters.")

    def _cached_flags(self, df) -> list:
        base = self._params()
        return ["✓" if os.path.exists(cache_path(
                    self.cache_dir, {**base, "contig": c, "start": int(a), "end": int(b)}))
                else "" for c, a, b in zip(df["contig"], df["start"], df["end"])]

    def _refresh_cached(self):
        df = self.seed_table.value
        if len(df):
            self.seed_table.patch({"cached": list(enumerate(self._cached_flags(df)))})

    def _on_seed(self, event):
        """Fill the seed fields from the picked row; recruit at once if it is cached."""
        if not event.new:
            return
        row = self.seed_table.value.iloc[event.new[0]]
        self.contig.value = row["contig"]
        self.start.value = int(row["start"])
        self.end.value = int(row["end"])
        if row["cached"]:
            self._on_recruit()
        else:
            self.status.object = (f"Seed **{row['region']}** (breakpoint {row['breakpoint']}) "
                                  "set. Press **Recruit** (minutes; cached afterwards).")

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
        if self.score_job is not None and not self.score_job.done:
            self.score_job.cancel()
        self.score = None
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
        self.score = None
        self._refresh_cached()
        self._rebuild_region()
        self._start_scoring(job.out_path)

    # ---- scoring (pair-HMM locus score) ---------------------------------------

    def _start_scoring(self, recruit_pkl: str):
        params = self._params()
        self.score_pane.object = self._score_html(running=True)
        self.score_job = ScoreJob(params, recruit_pkl, self.cache_dir).start()
        if self._score_poller is not None:
            self._score_poller.stop()
        self._score_poller = pn.state.add_periodic_callback(self._poll_score, period=500)

    def _poll_score(self):
        job = self.score_job
        if job is None:
            return
        self.score_log.value = job.log
        if not job.done:
            return
        if self._score_poller is not None:
            self._score_poller.stop()
            self._score_poller = None
        if job.error:
            self.score_pane.object = (f"<b>Scoring failed:</b> {job.error}. "
                                      "See the scoring log in section 3.")
            return
        if job.result is None:
            self.score_pane.object = "Scoring cancelled."
            return
        self.score = job.result
        self._render()

    def _score_html(self, running: bool = False) -> str:
        if running:
            return ("<div style='color:#666'>Scoring reads with the pair HMM… "
                    "(about 30 s per hundred reads; cached afterwards)</div>")
        ls = self.score
        if ls is None:
            return ""
        C = float(self.clamp.value)
        ds = ls.at(C)
        emis = ls.emission_at(C)
        up, down = ls.clamped_counts(C)
        fav = ls.asm1 if ds > 0 else ls.asm2 if ds < 0 else "neither"
        x = ls.exposure
        red = "#%02x%02x%02x" % RATIO_POS
        green = "#%02x%02x%02x" % RATIO_NEG
        return f"""
<div style="font-family:system-ui,sans-serif;font-size:13px;line-height:1.5;
            border:1px solid #ddd;border-radius:3px;padding:8px 12px;margin:4px 0">
  <div style="font-size:16px"><b>&Delta;S = {ds:+.2f} nats</b>
    &nbsp;&rarr; favours <b>{fav}</b>
    <span style="color:#888">&nbsp;(C = {C:g})</span></div>
  <div>emission <b>{emis:+.2f}</b>
       <span style="color:#888">({up} reads at +C, {down} at &minus;C)</span>
       &nbsp;&middot;&nbsp; exposure <b>{x.value:+.2f}</b>
       <span style="color:#888">&minus;2&lambda;(L<sub>eff</sub> {x.Leff1:,.0f} &minus; {x.Leff2:,.0f})</span></div>
  <div style="color:#666">{len(ls.reads)} reads scored &middot; {len(ls.dropped)} dropped at
       segment edges &middot; {len(ls.unscorable)} unscorable &middot; worst band edge
       {ls.worst_band_edge:.1e} &middot; &lambda; {x.lam:.3e} from {self.asm2.value}</div>
  <div style="display:flex;align-items:center;gap:8px;margin-top:6px;color:#666">
    <span>favours {ls.asm2}</span>
    <span style="flex:0 0 220px;height:9px;border-radius:2px;
                 background:linear-gradient(to right,{green},#d2d2ca,{red})"></span>
    <span>favours {ls.asm1}</span>
  </div>
</div>"""

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

    def _plot_box(self) -> dict:
        """Bounds for the dotplot so the whole report fits the browser window."""
        vw, vh = self.viewport.width, self.viewport.height
        if not vw or not vh:
            return {}
        margin = MARGIN_SIZE if vh >= SMALL_SCREEN else SMALL_MARGIN
        return dict(max_width=max(300, vw - CHROME_WIDTH - margin),
                    max_height=max(300, vh - CHROME_HEIGHT - margin),
                    margin_size=margin)

    def _render(self):
        rd = self.rd
        if rd is None:
            return
        ls = self.score
        ratios = {r.name: r.ratio for r in ls.reads} if ls is not None else None
        dropped = set(ls.dropped) if ls is not None else None
        self.plot.object = report_layout(
            rd, color_by=self.color_by.value, visible_layers=set(self.layers.value),
            ratios=ratios, dropped=dropped, C=float(self.clamp.value), **self._plot_box())
        if ls is not None:
            self.score_pane.object = self._score_html()
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

    @staticmethod
    def _pair(a, b):
        """Two widgets side by side inside the sidebar (default widths overflow it)."""
        a.width = b.width = 140
        return pn.Row(a, b)

    def sidebar(self):
        full = SIDEBAR_WIDTH - 40
        for widget in (self.base, self.label, self.contig, self.recruit_pid,
                       self.merge_gap, self.log, self.pid, self.plen, self.binsize,
                       self.layers, self.score_log, self.export_dir, self.export_binsize,
                       self.quast_run, self.seed_table,
                       self.export_log):
            widget.width = full
        return pn.Accordion(
            ("Seeds (QUAST misassemblies)", pn.Column(
                self.quast_run, self.seed_table,
                pn.Row(self.seed_reload, self.seed_note))),
            ("1 · Recruitment (minutes)", pn.Column(
                self.base, self.label, self._pair(self.asm1, self.asm2),
                self.contig, self._pair(self.start, self.end),
                self.recruit_pid, self.merge_gap,
                self._pair(self.ext_window, self.ext_max),
                self._pair(self.ext_rounds, self.max_rounds),
                pn.Row(self.recruit_btn, self.cancel_btn, self.busy),
                self.log)),
            ("2 · Fetch & classify (~1 s)", pn.Column(
                self.pid, self.plen,
                self._pair(self.best_pid_delta, self.pid_tie_eps),
                self.binsize, self.tie_warning)),
            ("3 · Layers & scoring log", pn.Column(
                pn.pane.Markdown("**layers**"), self.layers,
                self.score_log)),
            ("4 · Export (report_plot)", pn.Column(
                pn.pane.Markdown(
                    "Publication figures come from matplotlib `report_plot`, including "
                    "its k-mer dotplot — not from this page.", width=full),
                self.export_dir, self._pair(self.export_fmt, self.export_dpi),
                self.export_binsize, self.export_view_only,
                self.export_btn, self.export_log)),
            active=[0, 1], toggle=False, width=SIDEBAR_WIDTH,
        )

    def display_bar(self):
        """Colour mode and clamp, kept above the plot where they are always visible."""
        return pn.Row(pn.pane.Markdown("**colour by**", margin=(12, 5, 0, 10)),
                      self.color_by, self.clamp, width=640)

    def view(self):
        # Each pane scrolls on its own inside the window; the page itself never does.
        pane_css = {"height": "calc(100vh - 16px)", "overflow-y": "auto"}
        return pn.Row(
            pn.Column(self.sidebar(), width=SIDEBAR_WIDTH + 20, styles=pane_css),
            pn.Column(pn.Row(self.status, self.display_bar(), sizing_mode="stretch_width"),
                      self.score_pane, self.plot,
                      self.viewport, sizing_mode="stretch_width", styles=pane_css),
            sizing_mode="stretch_width")


def build_app(**overrides):
    return QaviewApp(**overrides).view()


if pn.state.served:
    build_app().servable(title="qaview")
