"""Background work for an open project: recruit every seed, merge into loci, score (§13.3).

Opening a project starts two queues:

- **recruitment** -- every seed not already on disk, `recruit_workers` at a time
  in worker processes (spawned, not forked: the server process has threads). A
  result is `(segments, reads)` pickled under `output/recruit/`, keyed by the
  seed, the data files and the recruitment parameters.
- **scoring** -- each distinct locus once. A finished seed's locus is
  `locus_id(segments, reads)`; seeds that recruit the same locus share a row.
  Scores live in `output/loci/`, the cache the batch notebooks write, so a locus
  the batch scored is loaded, not recomputed. Scoring runs in a subprocess
  (`qaview.score_job`), which forks `score_workers` processes of its own.

Nothing here is authoritative state: it is all rebuilt from those files, so a
restart loses only the jobs in flight.
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

from qaview.project import Project, Seed

# Bump when a recruiter change alters results, so older recruitments are redone.
RECRUIT_VERSION = 1


# ---- recruitment (worker side) ----------------------------------------------------

def recruit_key(project: Project, seed: Seed) -> str:
    payload = json.dumps({
        "version": RECRUIT_VERSION, "db": project.db,
        "asm": [[c.name, c.fasta, c.bam] for c in (project.asm1, project.asm2)],
        "seed": list(seed.interval), "params": project.recruitment,
    }, sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def recruit_one(job: dict) -> str:
    """Recruit one seed and pickle {"segments", "reads"} to job["out"]. Worker side."""
    import contextlib
    import io

    from qatool.SegmentStorage import Segment, SegmentCollection
    from qatool.SRRecruiter import Recruiter

    r = job["recruitment"]
    srr = Recruiter(
        db_path=job["db"], asm1=job["asm1"], asm2=job["asm2"],
        asm1_bam_path=job["bam1"], asm2_bam_path=job["bam2"],
        asm1_fasta_path=job["fa1"], asm2_fasta_path=job["fa2"],
        merge_gap=r["merge_gap"], pid_threshold=r["pid_threshold"],
        plen_threshold=r["plen_threshold"], best_pid_delta=r["best_pid_delta"],
        extension_window=r["extension_window"], ext_max=r["ext_max"],
        ext_rounds=r["ext_rounds"], max_rounds=r["max_rounds"],
        out_dir=None, verbose=False,
    )
    seed = SegmentCollection([Segment(*job["seed"])])
    with contextlib.redirect_stdout(io.StringIO()):     # the recruiter prints regardless
        segments, reads, _ = srr.run(seed)
    tmp = f"{job['out']}.{os.getpid()}.tmp"
    with open(tmp, "wb") as fh:
        pickle.dump({"segments": list(segments), "reads": set(reads)}, fh)
    os.replace(tmp, job["out"])
    return job["out"]


# ---- state --------------------------------------------------------------------------

@dataclass
class SeedState:
    seed: Seed
    path: str                   # recruitment pickle
    status: str = "queued"      # queued | done | failed
    locus: str | None = None
    error: str | None = None


@dataclass
class Locus:
    id: str
    recruit_path: str           # any one seed's pickle: all of them hold this locus
    seeds: list[str] = field(default_factory=list)
    n_reads: int = 0
    span: dict = field(default_factory=dict)       # assembly -> bp
    n_segments: dict = field(default_factory=dict)  # assembly -> count
    status: str = "queued"      # queued | scoring | scored | failed
    progress: float = 0.0       # scoring, 0..1
    score: object = None        # qamodel LocusScore
    error: str | None = None

    def row(self) -> dict:
        s = self.score
        row = {"id": self.id, "seeds": list(self.seeds), "n_seeds": len(self.seeds),
               "reads": self.n_reads, "span": self.span, "segments": self.n_segments,
               "status": self.status, "progress": round(self.progress, 3),
               "error": self.error, "delta_s": None}
        if s is not None:
            up, down = s.clamped_counts()
            row.update(delta_s=round(s.delta_s, 2), emission=round(s.emission, 2),
                       exposure=round(s.exposure.value, 2), scored=len(s.reads),
                       dropped=len(s.dropped), at_plus_C=up, at_minus_C=down, C=s.C)
        return row


class Pipeline:
    """Owns the queues for one project. Thread-safe; the server only reads snapshots."""

    def __init__(self, project: Project, python: str | None = None):
        self.project = project
        self.python = python or sys.executable
        os.makedirs(project.recruit_dir, exist_ok=True)
        os.makedirs(project.loci_dir, exist_ok=True)
        self.seeds = {s.name: SeedState(s, os.path.join(
            project.recruit_dir, f"recruit_{recruit_key(project, s)}.pkl"))
            for s in project.seeds}
        self.loci: dict[str, Locus] = {}
        self._lock = threading.RLock()
        self._wake = threading.Condition(self._lock)
        self._seq = 0
        self._changed: dict[str, int] = {}          # locus id -> seq of its last change
        self._front: list[str] = []                 # loci to score next, newest first
        self._pool: ProcessPoolExecutor | None = None
        self._score_proc: subprocess.Popen | None = None
        self._stopping = False
        self.started = time.time()

    # ---- lifecycle -------------------------------------------------------------------

    def start(self) -> "Pipeline":
        threading.Thread(target=self._recruit_all, daemon=True, name="recruit").start()
        threading.Thread(target=self._score_loop, daemon=True, name="score").start()
        return self

    def stop(self):
        with self._lock:
            self._stopping = True
            self._wake.notify_all()
        if self._pool is not None:
            # shutdown() only cancels queued seeds; running workers must be ended.
            procs = list(getattr(self._pool, "_processes", {}).values())
            self._pool.shutdown(wait=False, cancel_futures=True)
            for proc in procs:
                if proc.is_alive():
                    proc.terminate()
        proc = self._score_proc
        if proc is not None and proc.poll() is None:
            # The scoring job forks its own workers: end the whole process group.
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    # ---- recruitment -----------------------------------------------------------------

    def _job(self, st: SeedState) -> dict:
        p = self.project
        return {"db": p.db, "asm1": p.asm1.name, "asm2": p.asm2.name,
                "bam1": p.asm1.bam, "bam2": p.asm2.bam,
                "fa1": p.asm1.fasta, "fa2": p.asm2.fasta,
                "recruitment": p.recruitment, "seed": list(st.seed.interval),
                "out": st.path}

    def _recruit_all(self):
        # Seeds already on disk first, so a reopened project fills its table at once.
        todo = []
        for st in list(self.seeds.values()):
            if self._stopping:
                return
            if os.path.exists(st.path):
                self._seed_done(st)
            else:
                todo.append(st)
        if not todo:
            return
        import multiprocessing as mp
        self._pool = ProcessPoolExecutor(
            max_workers=max(1, min(self.project.resources["recruit_workers"], len(todo))),
            mp_context=mp.get_context("spawn"))
        for st in todo:
            fut = self._pool.submit(recruit_one, self._job(st))
            fut.add_done_callback(lambda f, st=st: self._on_recruited(st, f))

    def _on_recruited(self, st: SeedState, fut):
        if fut.cancelled():
            return
        exc = fut.exception()
        if exc is not None:
            with self._lock:
                st.status, st.error = "failed", f"{type(exc).__name__}: {exc}"
                self._bump(None)
            return
        self._seed_done(st)

    def _seed_done(self, st: SeedState):
        from qatool.plotters.likelihood import locus_id

        try:
            with open(st.path, "rb") as fh:
                rec = pickle.load(fh)
        except Exception as exc:                      # a torn or stale file: redo it
            with self._lock:
                st.status, st.error = "failed", f"unreadable recruitment: {exc}"
                self._bump(None)
            return
        segs, reads = rec["segments"], rec["reads"]
        lid = locus_id(segs, reads, self.project.asm1.name, self.project.asm2.name)
        with self._lock:
            st.status, st.locus = "done", lid
            loc = self.loci.get(lid)
            if loc is None:
                names = (self.project.asm1.name, self.project.asm2.name)
                loc = Locus(lid, st.path, n_reads=len(reads),
                            span={a: sum(s.end - s.start for s in segs if s.assembly == a)
                                  for a in names},
                            n_segments={a: sum(s.assembly == a for s in segs) for a in names})
                self.loci[lid] = loc
                self._wake.notify_all()
            loc.seeds.append(st.seed.name)
            self._bump(lid)

    # ---- scoring ---------------------------------------------------------------------

    def score_path(self, loc: Locus) -> str:
        from qatool.plotters.likelihood import score_cache_path

        p = self.project
        return score_cache_path(
            p.loci_dir, loc.id, p.asm1.name, p.asm2.name, db_path=p.db,
            bams={p.asm1.name: p.asm1.bam, p.asm2.name: p.asm2.bam},
            fastas={p.asm1.name: p.asm1.fasta, p.asm2.name: p.asm2.fasta},
            genome_fasta=p.genome_fasta, library_assembly=p.library_assembly)

    def prioritize(self, lid: str):
        """Score this locus next (the user opened it)."""
        with self._lock:
            if lid in self._front:
                self._front.remove(lid)
            self._front.insert(0, lid)
            self._wake.notify_all()

    def _next_to_score(self) -> Locus | None:
        for lid in self._front:
            loc = self.loci.get(lid)
            if loc is not None and loc.status == "queued":
                return loc
        return next((l for l in self.loci.values() if l.status == "queued"), None)

    def _score_loop(self):
        while True:
            with self._lock:
                loc = self._next_to_score()
                while loc is None and not self._stopping:
                    self._wake.wait(timeout=1.0)
                    loc = self._next_to_score()
                if self._stopping:
                    return
                loc.status = "scoring"
                self._bump(loc.id)
            try:
                self._score(loc)
            except Exception as exc:
                with self._lock:
                    loc.status, loc.error = "failed", f"{type(exc).__name__}: {exc}"
                    self._bump(loc.id)

    def _score(self, loc: Locus):
        path = self.score_path(loc)
        if not os.path.exists(path):
            self._run_score_job(loc)
            if self._stopping:
                return
        with open(path, "rb") as fh:
            score = pickle.load(fh)
        with self._lock:
            loc.score, loc.status, loc.progress = score, "scored", 1.0
            self._bump(loc.id)

    def _run_score_job(self, loc: Locus):
        env = dict(os.environ, PYTHONUNBUFFERED="1",
                   QAMODEL_SCORE_WORKERS=str(self.project.resources["score_workers"]))
        proc = subprocess.Popen(
            [self.python, "-m", "qaview.score_job", self.project.path, loc.recruit_path],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env,
            start_new_session=True)        # its own group, so stop() reaches its workers
        self._score_proc = proc
        tail = []
        for line in proc.stdout:
            line = line.rstrip()
            if line.startswith("scored "):
                done, total = line.split()[1].split("/")
                with self._lock:
                    loc.progress = int(done) / max(1, int(total))
                    self._bump(loc.id)
            else:
                tail = (tail + [line])[-20:]
        if proc.wait() != 0 and not self._stopping:
            raise RuntimeError(f"scoring exited with code {proc.returncode}: "
                               + " | ".join(tail[-3:]))

    # ---- snapshots for the server ----------------------------------------------------

    def _bump(self, lid: str | None):
        self._seq += 1
        if lid is not None:
            self._changed[lid] = self._seq

    def progress(self) -> dict:
        with self._lock:
            st = [s.status for s in self.seeds.values()]
            lo = [l.status for l in self.loci.values()]
            return {"seeds": len(st), "recruited": st.count("done"),
                    "recruit_failed": st.count("failed"), "loci": len(lo),
                    "scored": lo.count("scored"), "scoring": lo.count("scoring"),
                    "score_failed": lo.count("failed"), "seq": self._seq}

    def rows(self) -> list[dict]:
        with self._lock:
            return [l.row() for l in self.loci.values()]

    def changes_since(self, seq: int) -> tuple[int, list[dict]]:
        with self._lock:
            rows = [self.loci[lid].row() for lid, s in self._changed.items() if s > seq]
            return self._seq, rows

    def failed_seeds(self) -> list[dict]:
        with self._lock:
            return [{"seed": s.seed.name, "error": s.error}
                    for s in self.seeds.values() if s.status == "failed"]

    def recruitment(self, lid: str) -> dict:
        with self._lock:
            loc = self.loci[lid]
        with open(loc.recruit_path, "rb") as fh:
            return pickle.load(fh)
