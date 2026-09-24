"""The qaview web server for one open project (docs/qaview.md §13.5).

FastAPI serves a static page and a small JSON API; the page does the layout
(CSS grid) and asks for figures sized to its plot cell. Figures are Bokeh
`json_item`s from `glyphs.report_layout`, whose interaction is all client-side,
so there is no Bokeh server and no websocket. Background work belongs to
`pipeline.Pipeline`; this module only reads its snapshots.
"""

from __future__ import annotations

import asyncio
import json
import os
import threading
from collections import OrderedDict
from contextlib import asynccontextmanager

from fastapi import Body, FastAPI, HTTPException, Query
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

from qaview.glyphs import BAR_SIZE, MARGIN_SIZE, RATIO_NEG, RATIO_POS, report_layout
from qaview.pipeline import Pipeline
from qaview.project import Project
from qaview.region import (LAYER_COMMON, LAYER_OTHER, LAYER_SUPPLEMENTARY, LAYER_SUPPORT,
                           RegionData)

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

LAYERS = {"common": LAYER_COMMON, "unclassified": LAYER_OTHER,
          "support": LAYER_SUPPORT, "supplementary": LAYER_SUPPLEMENTARY}

# Pixels of the report that are not the dotplot: the y contig bar and right
# marginal across, the top marginal and x contig bar down, plus the readout and
# status lines above the grid and Bokeh's toolbar beside the dotplot.
PAD_W, PAD_H = 40, 64
SMALL_MARGIN, SMALL_CELL = 120, 700


def _hex(rgb) -> str:
    return "#%02x%02x%02x" % rgb


class RegionCache:
    """RegionData per (locus, fetch/classify settings): the ~1 s tier, kept warm."""

    def __init__(self, size: int = 8):
        self.size = size
        self._d: OrderedDict = OrderedDict()
        self._lock = threading.Lock()

    def get(self, key, build):
        with self._lock:
            if key in self._d:
                self._d.move_to_end(key)
                return self._d[key]
        rd = build()
        with self._lock:
            self._d[key] = rd
            while len(self._d) > self.size:
                self._d.popitem(last=False)
        return rd


class Stars:
    """Seeds the user starred, kept in the project's output directory.

    Keyed by seed name, not locus id: seed names come from the seed files and
    survive a change of recruitment parameters, which renames every locus. A
    locus row is starred when any of its seeds is.
    """

    def __init__(self, path: str):
        self.path = path
        self._lock = threading.Lock()
        try:
            with open(path) as fh:
                self.seeds = set(json.load(fh).get("starred", []))
        except (OSError, ValueError):
            self.seeds = set()

    def set(self, seeds, starred: bool) -> list:
        with self._lock:
            if starred:
                self.seeds |= set(seeds)
            else:
                self.seeds -= set(seeds)
            tmp = f"{self.path}.tmp"
            with open(tmp, "w") as fh:
                json.dump({"starred": sorted(self.seeds)}, fh, indent=1)
            os.replace(tmp, self.path)
            return sorted(self.seeds)


def create_app(project: Project, pipeline: Pipeline | None = None,
               start: bool = True) -> FastAPI:
    pipe = pipeline or Pipeline(project)
    regions = RegionCache()
    os.makedirs(project.output_dir, exist_ok=True)
    stars = Stars(os.path.join(project.output_dir, "qaview_stars.json"))

    @asynccontextmanager
    async def lifespan(_app):
        if start:
            pipe.start()
        yield
        pipe.stop()

    app = FastAPI(title="qaview", lifespan=lifespan)
    # Views are JSON of per-read columns, which compresses ~5-10x; over an ssh
    # tunnel that is most of the wait on a large locus.
    app.add_middleware(GZipMiddleware, minimum_size=1024, compresslevel=5)

    @app.middleware("http")
    async def revalidate_page(request, call_next):
        # Without this a browser keeps a stale app.js or style.css after an update
        # (heuristic caching) and the page silently runs old code. no-cache still
        # lets it reuse its copy -- it only asks first, and the ETag makes that cheap.
        response = await call_next(request)
        path = request.url.path
        if path == "/" or path.startswith("/static/"):
            response.headers["Cache-Control"] = "no-cache"
        return response
    app.state.pipeline = pipe
    app.mount("/static", StaticFiles(directory=STATIC), name="static")

    from bokeh.util.paths import bokehjs_path
    app.mount("/bokeh", StaticFiles(directory=os.path.join(bokehjs_path(), "js")),
              name="bokeh")

    @app.get("/")
    def index():
        return FileResponse(os.path.join(STATIC, "index.html"))

    @app.get("/api/project")
    def get_project():
        p = project
        return {
            "name": p.name, "path": p.path,
            "asm1": vars(p.asm1), "asm2": vars(p.asm2), "db": p.db,
            "library_assembly": p.library_assembly, "seed_sources": p.seed_sources,
            "recruitment": p.recruitment, "resources": p.resources,
            "output_dir": p.output_dir,
            "colors": {"ratio_pos": _hex(RATIO_POS), "ratio_neg": _hex(RATIO_NEG)},
            "layers": list(LAYERS),
        }

    @app.get("/api/loci")
    def get_loci():
        seq, rows = pipe.changes_since(-1)
        return {"seq": seq, "progress": pipe.progress(), "rows": pipe.rows(),
                "failed_seeds": pipe.failed_seeds()}

    @app.get("/api/events")
    async def events():
        """Server-sent events: queue progress and the rows that changed."""
        async def stream():
            seq = pipe.progress()["seq"]
            yield f"event: progress\ndata: {json.dumps(pipe.progress())}\n\n"
            while True:
                await asyncio.sleep(0.5)
                new, rows = pipe.changes_since(seq)
                if new == seq:
                    yield ": keep-alive\n\n"
                    continue
                seq = new
                payload = {"progress": pipe.progress(), "rows": rows,
                           "failed_seeds": pipe.failed_seeds()}
                yield f"event: update\ndata: {json.dumps(payload)}\n\n"
        return StreamingResponse(stream(), media_type="text/event-stream",
                                 headers={"Cache-Control": "no-cache"})

    @app.get("/api/stars")
    def get_stars():
        return {"seeds": sorted(stars.seeds)}

    @app.post("/api/stars")
    def set_stars(seeds: list[str] = Body(...), starred: bool = Body(...)):
        known = {s.name for s in project.seeds}
        unknown = [x for x in seeds if x not in known]
        if unknown:
            raise HTTPException(400, f"unknown seeds: {', '.join(unknown[:5])}")
        return {"seeds": stars.set(seeds, starred)}

    @app.post("/api/loci/{lid}/prioritize")
    def prioritize(lid: str):
        if lid not in pipe.loci:
            raise HTTPException(404, f"no locus {lid}")
        pipe.prioritize(lid)
        return {"ok": True}

    @app.get("/api/loci/{lid}/view")
    def view(lid: str, w: int = Query(900, ge=300), h: int = Query(700, ge=300),
             pid: float = 0.997, plen: float = 0.999,
             best_pid_delta: float = 0.001, pid_tie_eps: float = 0.0001,
             binsize: int = 0, color_by: str = "log-ratio", C: float = 200.0,
             layers: str = ",".join(LAYERS)):
        """The report for one locus, sized to a w x h cell, as a Bokeh json_item."""
        from bokeh.embed import json_item

        loc = pipe.loci.get(lid)
        if loc is None:
            raise HTTPException(404, f"no locus {lid}")
        if color_by not in ("class", "block", "log-ratio"):
            raise HTTPException(400, f"unknown colour mode {color_by!r}")
        p = project
        key = (lid, pid, plen, best_pid_delta, pid_tie_eps, binsize)

        def build():
            rec = pipe.recruitment(lid)
            return RegionData.build(
                rec["segments"], rec["reads"], p.db, p.asm1.bam, p.asm2.bam,
                asm1_fasta=p.asm1.fasta, asm2_fasta=p.asm2.fasta,
                asm1=p.asm1.name, asm2=p.asm2.name,
                pid_threshold=pid, plen_threshold=plen, binsize=binsize or None,
                best_pid_delta=best_pid_delta, pid_tie_eps=pid_tie_eps)

        rd = regions.get(key, build)
        score = loc.score
        ratios = {r.name: r.ratio for r in score.reads} if score is not None else None
        dropped = set(score.dropped) if score is not None else None
        margin = MARGIN_SIZE if min(w, h) >= SMALL_CELL else SMALL_MARGIN
        visible = {LAYERS[x] for x in layers.split(",") if x in LAYERS}
        layout = report_layout(
            rd, color_by=color_by, visible_layers=visible, ratios=ratios, dropped=dropped,
            C=C, max_width=w - BAR_SIZE - margin - PAD_W,
            max_height=h - BAR_SIZE - margin - PAD_H, margin_size=margin)
        out = {"item": json_item(layout), "locus": loc.row(),
               "tie_warning": pid_tie_eps >= best_pid_delta, "score": None}
        if score is not None:
            up, down = score.clamped_counts(C)
            x = score.exposure
            out["score"] = {
                "C": C, "delta_s": score.at(C), "emission": score.emission_at(C),
                "exposure": x.value, "at_plus_C": up, "at_minus_C": down,
                "scored": len(score.reads), "dropped": len(score.dropped),
                "unscorable": len(score.unscorable), "lam": x.lam,
                "Leff1": x.Leff1, "Leff2": x.Leff2,
                "worst_band_edge": score.worst_band_edge,
                "placements": getattr(score, "placements", 0),
                "pruned": getattr(score, "pruned", 0)}
        return JSONResponse(out)

    return app
