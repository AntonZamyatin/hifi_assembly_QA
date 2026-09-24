"""The pair-HMM locus score as a cached, cancellable background job.

Scoring sits between recruitment and display in cost: about 30 s for the HG002
seed locus's 106 reads on the compiled engine, longer for bigger units. So it
runs like recruitment -- in a subprocess, cached on disk -- and starts on its own
as soon as a recruitment result arrives.

The cache key is the recruitment key plus everything that changes a read's
log-ratio: lambda's reference, window sizing, band width, and the pair-HMM
parameters. The clamp C is deliberately *not* part of it. `LocusScore.at(C)`
re-evaluates Delta S from the stored per-read ratios, so C is a display setting
and moving its slider never rescores.
"""

import hashlib
import json
import os
import pickle
from dataclasses import asdict

from qaview.jobs import SubprocessJob
from qaview.recruit import recruit_key

WINDOW = dict(k=20, extra=50, width=100)


def score_config(params: dict) -> dict:
    """Everything besides the recruitment that changes a read's log-ratio.

    lambda is measured against asm2's primary reads over asm2's diploid length.
    That must be the same for every locus being compared, which holds as long as
    the assembly pair is fixed; it is reported in the panel so it is never implicit.
    """
    from qamodel.phmm import PhmmParams
    from qamodel.score import DEFAULT_PRUNE

    return dict(library_assembly=params["asm2"], genome_fasta=params["fa2"],
                phmm=asdict(PhmmParams()), prune=DEFAULT_PRUNE, **WINDOW)


def score_cache_path(cache_dir: str, params: dict) -> str:
    payload = json.dumps({"recruit": recruit_key(params), "score": score_config(params)},
                         sort_keys=True)
    return os.path.join(cache_dir,
                        f"score_{hashlib.sha1(payload.encode()).hexdigest()[:16]}.pkl")


def run_scoring(params: dict, recruit_pkl: str, out_path: str) -> None:
    """Score the recruited unit and pickle the LocusScore.  Child-side."""
    from qamodel import fasta_length, library_stats, score_locus

    with open(recruit_pkl, "rb") as fh:
        rec = pickle.load(fh)
    cfg = score_config(params)
    lib = library_stats(params["db"], cfg["library_assembly"],
                        fasta_length(cfg["genome_fasta"]))
    print(f"lambda {lib.lam:.4e} from {cfg['library_assembly']} "
          f"({lib.n_reads:,} reads, {lib.coverage:.1f}x per haplotype)")

    def progress(i, n):
        if i % 10 == 0 or i == n:
            print(f"scored {i}/{n} reads")

    result = score_locus(
        rec["reads"], rec["segments"], params["asm1"], params["asm2"],
        db_path=params["db"],
        bams={params["asm1"]: params["bam1"], params["asm2"]: params["bam2"]},
        fastas={params["asm1"]: params["fa1"], params["asm2"]: params["fa2"]},
        library=lib, progress=progress, **WINDOW)
    print(result.summary())
    with open(out_path, "wb") as fh:
        pickle.dump(result, fh)


class ScoreJob(SubprocessJob):
    """A running (or finished) scoring run.  Poll `lines`, `done`, `result`."""

    def __init__(self, params: dict, recruit_pkl: str, cache_dir: str,
                 python: str | None = None):
        self.params = params
        self.out_path = score_cache_path(cache_dir, params)
        self.result = None
        self.from_cache = False
        os.makedirs(cache_dir, exist_ok=True)
        self._param_file = self.out_path + ".params.json"
        super().__init__(["-m", "qaview.score_worker", self._param_file,
                          recruit_pkl, self.out_path], python)

    @property
    def label(self) -> str:
        return "scoring"

    def start(self) -> "ScoreJob":
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
