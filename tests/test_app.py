"""Structural tests for the qaview Panel app and recruitment job (Phase 3).

Synthetic fixtures only.  The behaviour worth pinning here is the *tier boundary*:
which widget triggers which cost.  Wiring a Tier-2 slider to `value` instead of
`value_throttled`, or letting a display toggle rebuild RegionData, would still work —
just badly — so it cannot be caught by eye.

Run: venv/bin/python tests/test_app.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

from bokeh.models import GridPlot

from fixtures import READS_A, READS_B, SEGS_A, SEGS_B, make_region_data, write_bam, write_db

from qaview.app import DEFAULTS, LAYER_LABELS, QaviewApp, dataset_paths
from qaview.recruit import RECRUIT_KEYS, RecruitJob, cache_path, recruit_key

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _app(tmp, **kw):
    return QaviewApp(cache_dir=os.path.join(tmp, "cache"), **kw)


def test_dataset_paths():
    p = dataset_paths("base", "lbl", "a1", "a2")
    check(p["db"] == "base/dbs/lbl.aln.sqlite", p["db"])
    check(p["bam1"] == "base/alignments/lbl.a1.sorted.bam", p["bam1"])
    check(p["fa2"] == "base/bgz/lbl.a2.fasta.bgz", p["fa2"])


def test_recruit_key_ignores_display_params(tmp):
    """Display settings must not invalidate the recruitment cache."""
    app = _app(tmp)
    base = app._params()
    k0 = recruit_key(base)

    check("pid_tie_eps" not in RECRUIT_KEYS, "pid_tie_eps is a display param")
    check("binsize" not in RECRUIT_KEYS, "binsize is a display param")

    app.pid.value = 0.95
    app.pid_tie_eps.value = 0.00001
    app.binsize.value = 4000
    check(recruit_key(app._params()) == k0,
          "changing display-only params must reuse the cached recruitment")

    app.merge_gap.value = 12345
    check(recruit_key(app._params()) != k0,
          "changing a recruitment param must produce a different cache key")


def test_cache_path_is_stable(tmp):
    app = _app(tmp)
    p = app._params()
    a, b = cache_path(tmp, p), cache_path(tmp, dict(p))
    check(a == b, "cache path must be deterministic")
    check(a.endswith(".pkl"), a)


def test_tier_wiring(tmp):
    """Tier 2 fires on release; Tier 3 fires immediately; neither re-recruits."""
    app = _app(tmp)

    def watches(widget, param, callback):
        """Is `callback` registered on this widget's `param`?

        Panel registers internal watchers on every parameter, so presence alone proves
        nothing — the check has to be for our own bound method.
        """
        for w in widget.param.watchers.get(param, {}).get("value", []):
            if getattr(w, "fn", None) == callback:
                return True
        return False

    for w, name in ((app.pid, "pid"), (app.plen, "plen")):
        check(watches(w, "value_throttled", app._on_tier2),
              f"Tier-2 slider {name} must rebuild on release (value_throttled)")
        check(not watches(w, "value", app._on_tier2),
              f"Tier-2 slider {name} must NOT rebuild on every drag step")

    for w, name in ((app.color_by, "color_by"), (app.layers, "layers")):
        check(watches(w, "value", app._on_tier3),
              f"Tier-3 control {name} must react immediately")
        check(not watches(w, "value", app._on_tier2),
              f"Tier-3 control {name} must not trigger a Tier-2 rebuild")


def test_tier3_does_not_rebuild_region(tmp):
    """A display change must reuse the cached RegionData, never re-read a BAM."""
    app = _app(tmp)
    app.rd = make_region_data(tmp, "app")
    calls = []
    app._rebuild_region = lambda: calls.append(1)     # would mean BAM/DB work
    app.layers.value = [0, 2]
    app.color_by.value = "block"
    check(not calls, "Tier-3 changes must not trigger a RegionData rebuild")
    check(app.plot.object is not None, "Tier-3 change must still redraw the figure")


def test_layer_toggle_hides_without_reindexing(tmp):
    """Hidden layers go to alpha 0; rows must stay in place for the CSR indexes."""
    app = _app(tmp)
    app.rd = make_region_data(tmp, "vis")
    n = len(app.rd.lines)

    app.layers.value = list(LAYER_LABELS)
    app._on_tier3()
    data = _main(app).renderers[0].data_source.data
    check(len(data["alpha"]) == n, "row count must not change with all layers on")
    all_on = sum(1 for a in data["alpha"] if a > 0)

    app.layers.value = []
    app._on_tier3()
    data = _main(app).renderers[0].data_source.data
    check(len(data["alpha"]) == n,
          "hiding layers must keep every row, or the hover CSR indexes misalign")
    check(all(a == 0 for a in data["alpha"]), "no layer selected -> nothing drawn")
    check(all_on > 0, "some lines should be visible when all layers are on")


def _main(app):
    grid = next(c for c in app.plot.object.children if isinstance(c, GridPlot))
    return {(r, c): ch for ch, r, c in grid.children}[(1, 1)]


def test_tie_warning(tmp):
    app = _app(tmp)
    app.best_pid_delta.value = 0.001

    app.pid_tie_eps.value = 0.0001
    app._check_tie()
    check(app.tie_warning.object == "", "a valid eps must not warn")

    app.pid_tie_eps.value = 0.001
    app._check_tie()
    check("pid_tie_eps" in (app.tie_warning.object or ""),
          "eps == delta makes asm1_best/asm2_best unreachable and must warn")

    app.pid_tie_eps.value = 0.01
    app._check_tie()
    check(app.tie_warning.object, "eps > delta must warn too")


def test_missing_inputs_reported(tmp):
    app = _app(tmp, base=os.path.join(tmp, "nope"))
    app._on_recruit()
    check("Missing inputs" in app.status.object,
          f"absent files must be reported, got {app.status.object!r}")
    check(app.job is None, "no job should start when inputs are missing")


def test_recruit_job_uses_cache(tmp):
    """A pre-existing cache file must short-circuit the subprocess entirely."""
    import pickle
    cache = os.path.join(tmp, "jobcache")
    os.makedirs(cache, exist_ok=True)
    params = {k: 1 for k in RECRUIT_KEYS}
    payload = {"segments": [], "reads": set()}
    with open(cache_path(cache, params), "wb") as fh:
        pickle.dump(payload, fh)

    job = RecruitJob(params, cache, python="/nonexistent/python").start()
    check(job.from_cache, "must load from cache")
    check(job.done and job.error is None, f"cached job should succeed, got {job.error}")
    check(job.result == payload, "cached payload must round-trip")


def test_recruit_worker_roundtrip(tmp):
    """The subprocess contract: params.json in, pickle out."""
    import json
    import pickle
    import subprocess

    bam1 = write_bam(os.path.join(tmp, "w1.bam"), READS_A)
    bam2 = write_bam(os.path.join(tmp, "w2.bam"), READS_B)
    rows = [{"assembly": "asm1", "read_name": "r1", "contig": "ctgA",
             "start": 2000, "end": 3500},
            {"assembly": "asm2", "read_name": "r1", "contig": "ctgB",
             "start": 1000, "end": 2500}]
    db = write_db(os.path.join(tmp, "w.sqlite"), rows)
    fa = os.path.join(tmp, "w.fa")
    with open(fa, "w") as fh:
        fh.write(">ctgA\n" + "ACGT" * 5000 + "\n>ctgB\n" + "ACGT" * 3000 + "\n")

    params = dict(db=db, bam1=bam1, bam2=bam2, fa1=fa, fa2=fa,
                  asm1="asm1", asm2="asm2", contig="ctgA", start=2000, end=6500,
                  pid_threshold=0.99, plen_threshold=0.5, best_pid_delta=0.001,
                  merge_gap=1000, extension_window=100, ext_max=10_000,
                  ext_rounds=1, max_rounds=2, out_dir=os.path.join(tmp, "scratch"))
    pf = os.path.join(tmp, "p.json")
    out = os.path.join(tmp, "o.pkl")
    with open(pf, "w") as fh:
        json.dump(params, fh)

    r = subprocess.run([sys.executable, "-m", "qaview.recruit_worker", pf, out],
                       capture_output=True, text=True, timeout=300)
    check(r.returncode == 0, f"worker failed: {r.stdout[-500:]} {r.stderr[-500:]}")
    if r.returncode == 0:
        with open(out, "rb") as fh:
            got = pickle.load(fh)
        check(set(got) == {"segments", "reads"}, f"unexpected payload keys {set(got)}")
        check(len(got["segments"]) > 0, "worker should recruit at least the seed")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_app_tests_")
    tests = [
        (test_dataset_paths, False),
        (test_recruit_key_ignores_display_params, True),
        (test_cache_path_is_stable, True),
        (test_tier_wiring, True),
        (test_tier3_does_not_rebuild_region, True),
        (test_layer_toggle_hides_without_reindexing, True),
        (test_tie_warning, True),
        (test_missing_inputs_reported, True),
        (test_recruit_job_uses_cache, True),
        (test_recruit_worker_roundtrip, True),
    ]
    for fn, needs_tmp in tests:
        print(f"{fn.__name__} ...")
        fn(tmp) if needs_tmp else fn()

    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
