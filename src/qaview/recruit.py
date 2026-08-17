"""Recruitment as a cached, cancellable background job.

Recruitment is the Tier-1 operation (minutes), so it must not block the server.  It
runs in a **subprocess** rather than a thread, which:

- sidesteps SQLite's thread affinity — `Recruiter.__init__` opens its connection at
  `SRRecruiter.py:45`, and here it is opened and used in the same process;
- lets progress be read from the child's stdout without hijacking `sys.stdout` in the
  server process;
- makes a hung or unwanted run killable.

Results are cached on disk keyed by a hash of the seed plus every parameter that
affects recruitment, so restarting the app or changing a *display* parameter never
re-recruits.
"""

import hashlib
import json
import os
import pickle

from qaview.jobs import SubprocessJob

# Parameters that change what recruitment produces.  Display-only settings are
# deliberately excluded so they never invalidate the cache.
RECRUIT_KEYS = (
    "db", "bam1", "bam2", "fa1", "fa2", "asm1", "asm2",
    "contig", "start", "end",
    "pid_threshold", "plen_threshold", "best_pid_delta",
    "merge_gap", "extension_window", "ext_max", "ext_rounds", "max_rounds",
)


def recruit_key(params: dict) -> str:
    payload = json.dumps({k: params[k] for k in RECRUIT_KEYS}, sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def cache_path(cache_dir: str, params: dict) -> str:
    return os.path.join(cache_dir, f"recruit_{recruit_key(params)}.pkl")


def run_recruitment(params: dict, out_path: str) -> None:
    """Run the recruiter and pickle (segments, reads) to out_path.  Child-side."""
    from qatool.SegmentStorage import Segment, SegmentCollection
    from qatool.SRRecruiter import Recruiter

    scratch = params.get("out_dir") or os.path.join(os.path.dirname(out_path), "recruit_scratch")
    os.makedirs(scratch, exist_ok=True)
    srr = Recruiter(
        db_path=params["db"], asm1=params["asm1"], asm2=params["asm2"],
        asm1_bam_path=params["bam1"], asm2_bam_path=params["bam2"],
        asm1_fasta_path=params["fa1"], asm2_fasta_path=params["fa2"],
        merge_gap=params["merge_gap"],
        pid_threshold=params["pid_threshold"], plen_threshold=params["plen_threshold"],
        best_pid_delta=params["best_pid_delta"],
        extension_window=params["extension_window"], ext_max=params["ext_max"],
        ext_rounds=params["ext_rounds"], max_rounds=params["max_rounds"],
        out_dir=scratch, verbose=True,
    )
    seed = SegmentCollection([Segment(params["asm1"], params["contig"],
                                      params["start"], params["end"])])
    segments, reads, _df = srr.run(seed)
    with open(out_path, "wb") as fh:
        pickle.dump({"segments": list(segments), "reads": set(reads)}, fh)


class RecruitJob(SubprocessJob):
    """A running (or finished) recruitment.  Poll `lines`, `done`, `result`."""

    def __init__(self, params: dict, cache_dir: str, python: str | None = None):
        self.params = params
        self.out_path = cache_path(cache_dir, params)
        self.result: dict | None = None
        self.from_cache = False
        os.makedirs(cache_dir, exist_ok=True)
        param_file = self.out_path + ".params.json"
        super().__init__(["-m", "qaview.recruit_worker", param_file, self.out_path], python)
        self._param_file = param_file

    @property
    def label(self) -> str:
        return "recruitment"

    def start(self) -> "RecruitJob":
        if os.path.exists(self.out_path):
            self.lines.append(f"[cache] reusing {os.path.basename(self.out_path)}")
            self._on_success()
            self.from_cache = True
            self.done = True
            return self
        with open(self._param_file, "w") as fh:
            json.dump(self.params, fh)
        super().start()
        return self

    def _on_success(self):
        with open(self.out_path, "rb") as fh:
            self.result = pickle.load(fh)
