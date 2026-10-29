"""The qaview web server (docs/qaview.md §13.5).

FastAPI serves a static page and a small JSON API; the page does the layout
(CSS grid) and asks for figures sized to its plot cell. Figures are Bokeh
`json_item`s from `glyphs.report_layout`, whose interaction is all client-side,
so there is no Bokeh server and no websocket. Background work belongs to
`pipeline.Pipeline`; this module only reads its snapshots.

The server may start with no project: the page then shows a form to open one
by its path on this machine. Everything that belongs to an open project lives
in a `Session`, replaced whole when another project is opened.
"""

from __future__ import annotations

import asyncio
import json
import os
import threading
from collections import OrderedDict
from contextlib import asynccontextmanager
from dataclasses import astuple, dataclass, fields, replace

from fastapi import Body, Depends, FastAPI, HTTPException, Query
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

from qaview.glyphs import (BACKGROUND_COLORS, BAR_BOX, MARGIN_SIZE, RATIO_NEG, RATIO_POS,
                           report_layout, round_colors, style_spec)
from qamodel.phmm.reads import FALLBACK_QUAL
from qaview.similarity import FLOOR as SIM_FLOOR
from qaview.pipeline import Pipeline
from qaview.project import (RECRUITMENT_CHOICES, RECRUITMENT_RANGES, Project, ProjectError,
                            load_project)
from qaview.region import (LAYER_COMMON, LAYER_OTHER, LAYER_PARTIAL, LAYER_SUPPLEMENTARY,
                           LAYER_SUPPORT, RegionData, arrange, read_alignments,
                           read_lengths)

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

LAYERS = {"common": LAYER_COMMON, "unclassified": LAYER_OTHER,
          "support": LAYER_SUPPORT, "supplementary": LAYER_SUPPLEMENTARY,
          "partial": LAYER_PARTIAL}

# "": the built layout; "x" / "y": that axis's segments along the other's diagonal
SORTS = ("", "x", "y")


@dataclass(frozen=True)
class RegionQuery:
    """The settings a locus's region is built for: query parameters of every
    endpoint that draws it or places things on it, with the page's defaults."""
    pid: float = 0.997
    plen: float = 0.999
    best_pid_delta: float = 0.001
    pid_tie_eps: float = 0.0001
    binsize: int = 0
    supp_plen: float = 0.3
    partial_pid: float = 0.99
    min_indel: int = 50
    sort: str = ""

# Pixels of the report that are not the dotplot: the y contig bar and right
# marginal across, the top marginal and x contig bar down, plus the readout and
# status lines above the grid and Bokeh's toolbar beside the dotplot.
PAD_W, PAD_H = 40, 64
SMALL_MARGIN, SMALL_CELL = 120, 700
TINY_MARGIN, TINY_CELL = 60, 450
MIN_CELL = 160              # the smallest cell a view is made for (static/app.js)


def _hex(rgb) -> str:
    return "#%02x%02x%02x" % rgb


class RegionCache:
    """RegionData per (locus, fetch/classify settings): the ~1 s tier, kept warm."""

    def __init__(self, size: int = 8):
        self.size = size
        self._d: OrderedDict = OrderedDict()
        self._lock = threading.Lock()
        self._building: dict = {}      # key -> lock held while it is built

    def get(self, key, build):
        with self._lock:
            if key in self._d:
                self._d.move_to_end(key)
                return self._d[key]
            building = self._building.setdefault(key, threading.Lock())
        # A second request for the same key waits for the first's build.
        with building:
            with self._lock:
                if key in self._d:
                    return self._d[key]
            rd = None
            try:
                rd = build()
            finally:
                with self._lock:
                    if rd is not None:
                        self._d[key] = rd
                        while len(self._d) > self.size:
                            self._d.popitem(last=False)
                    self._building.pop(key, None)
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


class Session:
    """An open project and what belongs to it: queues, region cache, stars."""

    def __init__(self, project: Project, pipeline: Pipeline | None = None):
        self.project = project
        self.pipe = pipeline or Pipeline(project)
        self.regions = RegionCache()
        os.makedirs(project.output_dir, exist_ok=True)
        self.stars = Stars(os.path.join(project.output_dir, "qaview_stars.json"))


RECENT_MAX = 30


def recent_path() -> str:
    base = os.environ.get("XDG_CONFIG_HOME") or os.path.expanduser("~/.config")
    return os.path.join(base, "qaview", "recent.json")


def recent_projects() -> list[dict]:
    try:
        with open(recent_path()) as fh:
            paths = json.load(fh)
    except (OSError, ValueError):
        return []
    return [{"path": p, "exists": os.path.exists(p)} for p in paths if isinstance(p, str)]


def remember_recent(path: str) -> None:
    paths = [path] + [r["path"] for r in recent_projects() if r["path"] != path]
    try:
        os.makedirs(os.path.dirname(recent_path()), exist_ok=True)
        with open(recent_path(), "w") as fh:
            json.dump(paths[:RECENT_MAX], fh, indent=1)
    except OSError:
        pass                                              # a convenience only


def create_app(project: Project | None = None, pipeline: Pipeline | None = None,
               start: bool = True) -> FastAPI:
    holder = {"session": Session(project, pipeline) if project is not None else None}
    switching = threading.Lock()

    def cur() -> Session:
        s = holder["session"]
        if s is None:
            raise HTTPException(409, "no project is open")
        return s

    def open_session(path: str) -> Session:
        """Load `path`, start its queues, and stop the previous project's."""
        new = Session(load_project(path))
        with switching:
            old, holder["session"] = holder["session"], new
        if start:
            new.pipe.start()
        if old is not None:
            old.pipe.stop()
        remember_recent(new.project.path)
        return new

    @asynccontextmanager
    async def lifespan(_app):
        s = holder["session"]
        if s is not None:
            if start:
                s.pipe.start()
            remember_recent(s.project.path)
        yield
        if holder["session"] is not None:
            holder["session"].pipe.stop()

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
    app.state.holder = holder
    app.mount("/static", StaticFiles(directory=STATIC), name="static")

    from bokeh.util.paths import bokehjs_path
    app.mount("/bokeh", StaticFiles(directory=os.path.join(bokehjs_path(), "js")),
              name="bokeh")

    @app.get("/")
    def index():
        return FileResponse(os.path.join(STATIC, "index.html"))

    @app.get("/api/project")
    def get_project():
        s = holder["session"]
        if s is None:
            return {"open": False, "recent": recent_projects()}
        p = s.project
        return {
            "open": True, "recent": recent_projects(),
            "name": p.name, "path": p.path,
            "asm1": vars(p.asm1), "asm2": vars(p.asm2), "db": p.db,
            "library_assembly": p.library_assembly, "seed_sources": p.seed_sources,
            "recruitment": p.recruitment, "resources": p.resources,
            "recruitment_ranges": {k: [t.__name__, lo, hi]
                                   for k, (t, lo, hi) in RECRUITMENT_RANGES.items()},
            "recruitment_choices": {k: list(v) for k, v in RECRUITMENT_CHOICES.items()},
            "output_dir": p.output_dir,
            "colors": {"ratio_pos": _hex(RATIO_POS), "ratio_neg": _hex(RATIO_NEG)},
            "background_colors": list(BACKGROUND_COLORS),
            "fallback_qual": FALLBACK_QUAL,
            "layers": list(LAYERS),
        }

    @app.post("/api/open")
    def open_project(path: str = Body(..., embed=True)):
        """Open the project file at `path` (on this machine), replacing the current one."""
        path = os.path.abspath(os.path.expanduser(path.strip()))
        if not path.endswith(".toml"):
            raise HTTPException(400, f"not a project file (.toml): {path}")
        if not os.path.exists(path):
            raise HTTPException(400, f"no such file: {path}")
        try:
            s = open_session(path)
        except ProjectError as exc:
            raise HTTPException(400, str(exc)) from exc
        return {"ok": True, "name": s.project.name, "path": s.project.path}

    @app.get("/api/loci")
    def get_loci():
        pipe = cur().pipe
        progress = pipe.progress()          # its seq is read before the rows
        return {"seq": progress["seq"], "progress": progress, "rows": pipe.rows(),
                "failed_seeds": pipe.failed_seeds()}

    @app.get("/api/events")
    async def events():
        """Server-sent events: queue progress and the rows that changed.

        A stream follows the project open when it started; when another is
        opened it ends, and the page reloads onto the new one.
        """
        pipe = cur().pipe

        async def stream():
            seq = pipe.progress()["seq"]
            yield f"event: progress\ndata: {json.dumps(pipe.progress())}\n\n"
            while True:
                await asyncio.sleep(0.5)
                if holder["session"] is None or holder["session"].pipe is not pipe:
                    return
                new, rows = pipe.changes_since(seq)
                if new == seq:
                    yield ": keep-alive\n\n"
                    continue
                seq = new
                payload = {"seq": new, "progress": pipe.progress(), "rows": rows,
                           "failed_seeds": pipe.failed_seeds()}
                yield f"event: update\ndata: {json.dumps(payload)}\n\n"
        return StreamingResponse(stream(), media_type="text/event-stream",
                                 headers={"Cache-Control": "no-cache"})

    @app.get("/api/stars")
    def get_stars():
        return {"seeds": sorted(cur().stars.seeds)}

    @app.post("/api/stars")
    def set_stars(seeds: list[str] = Body(...), starred: bool = Body(...)):
        s = cur()
        known = {x.name for x in s.project.seeds}
        unknown = [x for x in seeds if x not in known]
        if unknown:
            raise HTTPException(400, f"unknown seeds: {', '.join(unknown[:5])}")
        return {"seeds": s.stars.set(seeds, starred)}

    # ---- remembered alignments (qaview.remember) ------------------------------------

    def _layout(s: Session, lid: str, q: RegionQuery | None = None):
        """(segs_info1, segs_info2, flips1, flips2) of the locus as drawn. The built
        layout needs no reads; a sorted one is the view's own region (`q`: the
        view's query, cached by the view request that drew it)."""
        from qatool.plotters.axis import _build_axis_map

        if q and q.sort:
            rd = region(s, lid, q)
            return (rd.asm1_segs_info, rd.asm2_segs_info,
                    rd.flips_of(rd.asm1), rd.flips_of(rd.asm2))
        segs = sorted(s.pipe.recruitment(lid)["segments"], key=lambda g: (g.contig, g.start, g.end))
        return (_build_axis_map(segs, s.project.asm1.name)[0],
                _build_axis_map(segs, s.project.asm2.name)[0], None, None)

    @app.get("/api/loci/{lid}/remembered")
    def remembered_here(lid: str, q: RegionQuery = Depends()):
        """Remembered alignments whose stretches lie in this locus, as plot rectangles
        in the view drawn with these settings."""
        from qaview import remember

        s = cur()
        if lid not in s.pipe.loci:
            raise HTTPException(404, f"no locus {lid}")
        s1, s2, f1, f2 = _layout(s, lid, q)
        return remember.place(remember.all_meta(s.project.output_dir), s1, s2,
                              s.project.asm1.name, s.project.asm2.name, f1, f2)

    @app.post("/api/remembered")
    def remember_this(lid: str = Body(...), alignment: dict = Body(...),
                      note: str = Body(""), settings: dict = Body({}),
                      view: dict = Body({})):
        from qaview import remember

        s = cur()
        loc = s.pipe.loci.get(lid)
        if loc is None:
            raise HTTPException(404, f"no locus {lid}")
        try:
            meta = remember.save(s.project.output_dir, alignment, note=note,
                                 seeds=list(loc.seeds), locus=lid, settings=settings)
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc
        try:
            # typed as the view's query was, so the key finds the view's cached region
            args = (RegionQuery(**{f.name: type(f.default)(view[f.name]) for f in fields(RegionQuery)})
                    if view.get("sort") else None)
        except (KeyError, TypeError, ValueError):
            args = None
        s1, s2, f1, f2 = _layout(s, lid, args)
        placed = remember.place([meta], s1, s2, s.project.asm1.name, s.project.asm2.name, f1, f2)
        return {**meta, "rect": placed[0] if placed else None}

    @app.get("/api/remembered/{ident}")
    def remembered_one(ident: str):
        from qaview import remember

        try:
            return remember.load(cur().project.output_dir, ident)
        except FileNotFoundError as exc:
            raise HTTPException(404, f"no remembered alignment {ident}") from exc

    @app.delete("/api/remembered/{ident}")
    def forget_one(ident: str):
        from qaview import remember

        try:
            remember.forget(cur().project.output_dir, ident)
        except FileNotFoundError as exc:
            raise HTTPException(404, f"no remembered alignment {ident}") from exc
        return {"ok": True}

    # ---- adding a seed ------------------------------------------------------------------

    @app.post("/api/seeds")
    def add_seed(assembly: str = Body(...), region: str = Body(...), name: str = Body("")):
        """A new seed on `assembly` at `region` (contig:start-end, 1-based inclusive):
        saved to the output folder's qaview_seeds.json and recruited at once."""
        from qaview.project import new_seed, save_added_seed

        sess = cur()
        with switching:                    # one check-and-add at a time
            try:
                seed = new_seed(sess.project, assembly, region, name)
            except ValueError as exc:
                raise HTTPException(400, str(exc)) from exc
            save_added_seed(sess.project.output_dir, seed)
            sess.pipe.add_seed(seed)
        return {"name": seed.name, "assembly": seed.assembly, "contig": seed.contig,
                "start": seed.start, "end": seed.end}

    # ---- re-recruitment ----------------------------------------------------------------

    @app.get("/api/loci/{lid}/seeds")
    def locus_seeds(lid: str):
        """The locus's seeds, each with the recruitment parameters it was recruited with."""
        s = cur()
        try:
            seeds = s.pipe.seed_info(lid)
        except KeyError as exc:             # never was, or re-recruited away just now
            raise HTTPException(404, f"no locus {lid}") from exc
        return {"defaults": s.project.recruitment, "seeds": seeds}

    @app.post("/api/recruit")
    def rerecruit(seeds: list[str] = Body(...), params: dict = Body(...)):
        """Recruit these seeds again with `params`; the result arrives as row updates."""
        if not seeds:
            raise HTTPException(400, "no seeds given")
        try:
            return cur().pipe.rerecruit(seeds, params)
        except KeyError as exc:
            raise HTTPException(400, f"unknown seeds: {exc.args[0]}") from exc
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc

    @app.post("/api/loci/{lid}/prioritize")
    def prioritize(lid: str):
        pipe = cur().pipe
        if lid not in pipe.loci:
            raise HTTPException(404, f"no locus {lid}")
        pipe.prioritize(lid)
        return {"ok": True}

    def region(s: Session, lid: str, q: RegionQuery) -> RegionData:
        """The locus reduced to arrays for these fetch/classify settings (cached).
        `q.sort` "x" or "y" lays that axis's segments along the other's diagonal,
        flipping those the reads run backwards on (`region.arrange`)."""
        p = s.project
        if q.sort not in SORTS:
            raise HTTPException(400, f"unknown sort {q.sort!r}")
        key = (lid, *astuple(q))
        if q.sort:
            return s.regions.get(key, lambda: arrange(region(s, lid, replace(q, sort="")), q.sort))

        def build():
            rec = s.pipe.recruitment(lid)
            return RegionData.build(
                rec["segments"], rec["reads"], p.db, p.asm1.bam, p.asm2.bam,
                asm1_fasta=p.asm1.fasta, asm2_fasta=p.asm2.fasta,
                asm1=p.asm1.name, asm2=p.asm2.name,
                pid_threshold=q.pid, plen_threshold=q.plen, binsize=q.binsize or None,
                best_pid_delta=q.best_pid_delta, pid_tie_eps=q.pid_tie_eps,
                supp_plen_threshold=q.supp_plen, partial_pid_threshold=q.partial_pid,
                min_indel=q.min_indel)

        return s.regions.get(key, build)

    # ---- read browser -----------------------------------------------------------------

    @app.post("/api/reads/lengths")
    def reads_lengths(names: list[str] = Body(..., embed=True)):
        """Read -> length, for the reads of a clicked marginal bin."""
        p = cur().project
        return read_lengths(p.db, names, [p.asm1.name, p.asm2.name])

    @app.get("/api/reads/alignments")
    def reads_alignments(name: str, pid: float = 0.997, plen: float = 0.999,
                         best_pid_delta: float = 0.001):
        """Every alignment of one read in both candidates, marked (region.read_alignments)."""
        p = cur().project
        return read_alignments(p.db, name, [p.asm1.name, p.asm2.name], pid, plen, best_pid_delta)

    @app.get("/api/loci/{lid}/align")
    def align(lid: str, x0: float, x1: float, y0: float, y1: float,
              q: RegionQuery = Depends(), reads: bool = False, color_by: str = "log-ratio"):
        """Align the stretches under a rectangle on the dotplot (docs/qaview.md §13.4).

        With `reads`, also the reads on each stretch in the alignment's columns,
        coloured and banded as the marginals are for `color_by`.
        """
        import pysam

        from qaview.align import align_rect
        from qaview.glyphs import class_stacking, ratio_stacking, round_stacking
        from qaview.pileup import pileup

        s = cur()
        loc = s.pipe.loci.get(lid)
        if loc is None:
            raise HTTPException(404, f"no locus {lid}")
        rd = region(s, lid, q)
        p = s.project
        with pysam.FastaFile(p.asm1.fasta) as f1, pysam.FastaFile(p.asm2.fasta) as f2:
            try:
                out = align_rect(rd, f1, f2, x0, x1, y0, y1)
            except ValueError as exc:
                raise HTTPException(400, str(exc)) from exc
        if reads:
            score = loc.score
            rounds = s.pipe.recruitment(lid).get("rounds") if color_by == "round" else None
            stacking = (ratio_stacking(rd, {r.name: r.ratio for r in score.reads})
                        if color_by == "log-ratio" and score is not None else
                        round_stacking(rd, rounds) if rounds else class_stacking(rd))
            out["reads"] = pileup(
                rd, out, {p.asm1.name: p.asm1.bam, p.asm2.name: p.asm2.bam},
                {p.asm1.name: p.asm1.fasta, p.asm2.name: p.asm2.fasta}, p.db, stacking,
                pid=q.pid, plen=q.plen, best_pid_delta=q.best_pid_delta)
        return out

    @app.get("/api/loci/{lid}/view")
    def view(lid: str, w: int = Query(900, ge=MIN_CELL), h: int = Query(700, ge=MIN_CELL),
             q: RegionQuery = Depends(), color_by: str = "log-ratio", C: float = 200.0,
             layers: str = ",".join(LAYERS), reads: str = "hover", background: bool = True):
        """The report for one locus, sized to a w x h cell, as a Bokeh json_item.

        `reads`: "all" draws every read line, "hover" only the reads of a hovered
        marginal bin. `background`: the segments' k-mer similarity underneath, from
        `pid` up, once computed (it is started here; the row's `background` says
        when it is ready).
        """
        from bokeh.embed import json_item

        s = cur()
        loc = s.pipe.loci.get(lid)
        if loc is None:
            raise HTTPException(404, f"no locus {lid}")
        if color_by not in ("class", "block", "log-ratio", "round"):
            raise HTTPException(400, f"unknown colour mode {color_by!r}")
        if reads not in ("all", "hover"):
            raise HTTPException(400, f"unknown reads mode {reads!r}")
        bg = None
        if background:
            s.pipe.request_background(lid)
            bg = s.pipe.background(lid)
        # The background shows what the recruitment followed: drawn from the
        # locus's own recruitment pid, whatever the display pid (lowering that to
        # see more reads must not flood the background). The cache keeps cells
        # from similarity.FLOOR up, so a lower recruitment pid draws from there.
        bg_floor = max((loc.params or s.project.recruitment)["pid_threshold"], SIM_FLOOR)
        rd = region(s, lid, q)
        rounds = s.pipe.recruitment(lid).get("rounds")    # None: recruited before rounds were kept
        score = loc.score
        ratios = {r.name: r.ratio for r in score.reads} if score is not None else None
        dropped = set(score.dropped) if score is not None else None
        small = min(w, h)
        margin = (MARGIN_SIZE if small >= SMALL_CELL else
                  SMALL_MARGIN if small >= TINY_CELL else TINY_MARGIN)
        visible = {LAYERS[x] for x in layers.split(",") if x in LAYERS}
        layout = report_layout(
            rd, color_by=color_by, visible_layers=visible, ratios=ratios, dropped=dropped,
            C=C, max_width=w - BAR_BOX - margin - PAD_W,
            max_height=h - BAR_BOX - margin - PAD_H, margin_size=margin,
            background=bg, floor=bg_floor, reads=reads, rounds=rounds)
        out = {"item": json_item(layout), "locus": loc.row(),
               "background": {"status": loc.background, "error": loc.background_error,
                              "drawn": bg is not None, "floor": bg_floor,
                              "cells": int((bg["identity"] >= bg_floor).sum()) if bg is not None else 0,
                              "grid": [bg["n1"], bg["n2"]] if bg is not None else None,
                              "binsize": bg["binsize"] if bg is not None else None,
                              "sampling": bg["settings"]["sampling_scale"] if bg is not None else None},
               "tie_warning": q.pid_tie_eps >= q.best_pid_delta, "score": None,
               "rounds": round_colors(max(rounds.values(), default=1)) if rounds else None,
               # what the page needs to restyle without coming back: colouring,
               # C, layers, reads mode and the background are all applied there
               "style": {**style_spec(), "layers": dict(LAYERS)}}
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
                "pruned": getattr(score, "pruned", 0),
                "fixed_quality": len(getattr(score, "fixed_quality", ())),
                # per read, for Delta S at another C on the page (+/-inf as 1e308:
                # JSON has no infinity, and any C clamps it the same)
                "ratios": [max(-1e308, min(1e308, r.ratio)) for r in score.reads]}
        return JSONResponse(out)

    return app
