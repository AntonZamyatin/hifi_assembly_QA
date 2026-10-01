"""The project file: one comparison of two assembly candidates (docs/qaview.md §13.2).

A TOML file names the two candidates (FASTA and BAM each), the alignment
database, where the seeds come from and where results go. Paths are relative to
the file. Everything the app shows is derived from it, so opening the file is the
whole of the user's setup.

Seeds are source-neutral: each `[[seeds]]` entry is read by the loader for its
`format` into `Seed` records, and a new source is a new loader in `LOADERS`.
"""

from __future__ import annotations

import csv
import os
import tomllib
from collections import defaultdict
from dataclasses import dataclass, field

# The notebooks' values (examples/scripts/generate_notebooks_hg002_hic.py), so a
# project without [recruitment] reproduces the batch.
RECRUITMENT_DEFAULTS = dict(
    pid_threshold=0.997, plen_threshold=0.999, best_pid_delta=0.001,
    merge_gap=50_000, extension_window=100, ext_max=300_000,
    ext_rounds=100, max_rounds=100, extension="legacy",
)
RESOURCE_DEFAULTS = dict(recruit_workers=32, score_workers=32)
# What each recruitment parameter may be: (type, lowest, highest).
RECRUITMENT_RANGES = dict(
    pid_threshold=(float, 0.0, 1.0), plen_threshold=(float, 0.0, 1.0),
    best_pid_delta=(float, 0.0, 1.0), merge_gap=(int, 0, 10_000_000),
    extension_window=(int, 1, 100_000), ext_max=(int, 0, 100_000_000),
    ext_rounds=(int, 0, 10_000), max_rounds=(int, 1, 10_000),
)
# Parameters that take one of a few names (qatool.SRRecruiter.EXTENSIONS).
RECRUITMENT_CHOICES = dict(extension=("legacy", "anchor"))


class ProjectError(ValueError):
    """The project file is malformed or names files that do not exist."""


def validate_recruitment(params: dict) -> dict:
    """Recruitment parameters checked and coerced to their types; ValueError if not.

    Only the keys given are returned, so this checks partial overrides as well.
    """
    out, problems = {}, []
    for k, v in params.items():
        if k in RECRUITMENT_CHOICES:
            if v not in RECRUITMENT_CHOICES[k]:
                problems.append(f"{k} must be one of {', '.join(RECRUITMENT_CHOICES[k])}: {v!r}")
            else:
                out[k] = v
            continue
        if k not in RECRUITMENT_RANGES:
            problems.append(f"unknown parameter {k}")
            continue
        kind, lo, hi = RECRUITMENT_RANGES[k]
        try:
            x = float(v)
            if kind is int:
                if x != int(x):
                    raise ValueError
                x = int(x)
        except (TypeError, ValueError):
            problems.append(f"{k} must be {'a whole number' if kind is int else 'a number'}: {v!r}")
            continue
        if not lo <= x <= hi:
            problems.append(f"{k} must be between {lo:,} and {hi:,}: {v!r}")
            continue
        out[k] = x
    if problems:
        raise ValueError("; ".join(problems))
    return out


@dataclass(frozen=True)
class Candidate:
    name: str
    fasta: str
    bam: str


@dataclass(frozen=True)
class Seed:
    """One recruitment seed: an interval on one candidate, from one source."""

    name: str                   # unique within the project
    source: str                 # the [[seeds]] entry's format, e.g. "quast"
    assembly: str
    contig: str
    start: int                  # 0-based half-open
    end: int
    meta: dict = field(default_factory=dict, compare=False, hash=False)

    @property
    def interval(self) -> tuple:
        return (self.assembly, self.contig, self.start, self.end)


@dataclass
class Project:
    path: str
    name: str
    asm1: Candidate
    asm2: Candidate
    db: str
    library_assembly: str
    seeds: list[Seed]
    recruitment: dict
    resources: dict
    output_dir: str
    seed_sources: list[dict]

    @property
    def recruit_dir(self) -> str:
        return os.path.join(self.output_dir, "recruit")

    @property
    def loci_dir(self) -> str:
        return os.path.join(self.output_dir, "loci")

    def candidate(self, name: str) -> Candidate:
        for c in (self.asm1, self.asm2):
            if c.name == name:
                return c
        raise KeyError(name)

    @property
    def genome_fasta(self) -> str:
        """lambda's genome: the library assembly's FASTA."""
        return self.candidate(self.library_assembly).fasta


# ---- seed loaders ---------------------------------------------------------------

def load_quast(path: str, assembly: str) -> list[Seed]:
    """misassemblies.tsv from examples/quast/extract_misassemblies.py.

    Names follow the notebooks, `<contig>_<type>_<i:03d>` with `i` counting rows
    of one (contig, type) in file order, so a seed matches its notebook figure
    and batch score.
    """
    out, seen = [], defaultdict(int)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seen[(r["contig"], r["type"])] += 1
            i = seen[(r["contig"], r["type"])]
            bp0, bp1 = int(r["bp_start"]), int(r["bp_end"])
            meta = {"type": r["type"],
                    "breakpoint": f"{bp0:,}" if bp1 - bp0 <= 1 else f"{bp0:,}-{bp1:,}"}
            if r.get("inconsistency_bp"):
                meta["inconsistency_bp"] = int(r["inconsistency_bp"])
            out.append(Seed(f"{r['contig']}_{r['type']}_{i:03d}", "quast", assembly,
                            r["contig"], int(r["start"]), int(r["end"]), meta))
    return out


def load_bed(path: str, assembly: str) -> list[Seed]:
    """contig, start, end[, name, ...]: tab-separated, 0-based half-open.

    Blank lines and lines starting with '#', 'track' or 'browser' are skipped.
    Without a name column a seed is called `<contig>_<start>_<end>`.
    """
    out = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                raise ProjectError(f"{path}:{n}: expected contig, start, end")
            contig, start, end = f[0], int(f[1]), int(f[2])
            if end <= start:
                raise ProjectError(f"{path}:{n}: empty interval {start}-{end}")
            name = f[3] if len(f) > 3 and f[3] else f"{contig}_{start}_{end}"
            out.append(Seed(name, "bed", assembly, contig, start, end))
    return out


LOADERS = {"quast": load_quast, "bed": load_bed}


# ---- loading ----------------------------------------------------------------------

def _need(table: dict, key: str, where: str):
    if key not in table:
        raise ProjectError(f"{where}: missing '{key}'")
    return table[key]


def load_project(path: str, check_files: bool = True) -> Project:
    """Parse and validate a project file. Raises ProjectError with every problem found."""
    path = os.path.abspath(path)
    base = os.path.dirname(path)
    with open(path, "rb") as fh:
        try:
            t = tomllib.load(fh)
        except tomllib.TOMLDecodeError as exc:
            raise ProjectError(f"{path}: {exc}") from exc

    def p(rel: str) -> str:
        return os.path.normpath(os.path.join(base, os.path.expanduser(rel)))

    cands = []
    for key in ("asm1", "asm2"):
        c = _need(t, key, path)
        cands.append(Candidate(_need(c, "name", f"[{key}]"),
                               p(_need(c, "fasta", f"[{key}]")),
                               p(_need(c, "bam", f"[{key}]"))))
    asm1, asm2 = cands
    if asm1.name == asm2.name:
        raise ProjectError("asm1 and asm2 must have different names")

    data = _need(t, "data", path)
    db = p(_need(data, "db", "[data]"))
    library = data.get("library_assembly", asm2.name)
    if library not in (asm1.name, asm2.name):
        raise ProjectError(f"[data] library_assembly '{library}' is neither candidate")

    unknown = set(t.get("recruitment", {})) - set(RECRUITMENT_DEFAULTS)
    if unknown:
        raise ProjectError(f"[recruitment] unknown keys: {', '.join(sorted(unknown))}")
    recruitment = {**RECRUITMENT_DEFAULTS, **t.get("recruitment", {})}
    for k, names in RECRUITMENT_CHOICES.items():
        if recruitment[k] not in names:
            raise ProjectError(f"[recruitment] {k} must be one of {', '.join(names)}: "
                               f"{recruitment[k]!r}")
    unknown = set(t.get("resources", {})) - set(RESOURCE_DEFAULTS)
    if unknown:
        raise ProjectError(f"[resources] unknown keys: {', '.join(sorted(unknown))}")
    resources = {**RESOURCE_DEFAULTS, **t.get("resources", {})}

    output_dir = p(_need(_need(t, "output", path), "dir", "[output]"))

    sources = t.get("seeds", [])
    if isinstance(sources, dict):
        sources = [sources]
    if not sources:
        raise ProjectError("no [[seeds]] sources")
    problems = []
    for key, c in (("asm1", asm1), ("asm2", asm2)):
        for what in ("fasta", "bam"):
            f = getattr(c, what)
            if check_files and not os.path.exists(f):
                problems.append(f"[{key}] {what} not found: {f}")
    if check_files and not os.path.exists(db):
        problems.append(f"[data] db not found: {db}")

    seeds, source_info = [], []
    for i, src in enumerate(sources, 1):
        where = f"[[seeds]] #{i}"
        fmt = _need(src, "format", where)
        if fmt not in LOADERS:
            raise ProjectError(f"{where}: unknown format '{fmt}' "
                               f"(known: {', '.join(sorted(LOADERS))})")
        spath = p(_need(src, "path", where))
        assembly = src.get("assembly", asm1.name)
        if assembly not in (asm1.name, asm2.name):
            raise ProjectError(f"{where}: assembly '{assembly}' is neither candidate")
        if not os.path.exists(spath):
            problems.append(f"{where}: seed file not found: {spath}")
            continue
        got = LOADERS[fmt](spath, assembly)
        seeds += got
        source_info.append({"format": fmt, "path": spath, "assembly": assembly,
                            "seeds": len(got)})
    if problems:
        raise ProjectError("\n".join(problems))

    return Project(path=path, name=t.get("name", os.path.splitext(os.path.basename(path))[0]),
                   asm1=asm1, asm2=asm2, db=db, library_assembly=library,
                   seeds=_dedupe(seeds), recruitment=recruitment, resources=resources,
                   output_dir=output_dir, seed_sources=source_info)


def _dedupe(seeds: list[Seed]) -> list[Seed]:
    """One seed per interval; names from several sources are kept in meta['aliases'].

    A name used twice for different intervals gets a numeric suffix, so names stay
    unique -- they key the recruitment results and the table.
    """
    by_interval: dict[tuple, Seed] = {}
    aliases: dict[tuple, list[str]] = defaultdict(list)
    for s in seeds:
        if s.interval in by_interval:
            aliases[s.interval].append(s.name)
        else:
            by_interval[s.interval] = s
    out, used = [], defaultdict(int)
    for iv, s in by_interval.items():
        used[s.name] += 1
        name = s.name if used[s.name] == 1 else f"{s.name}.{used[s.name]}"
        meta = dict(s.meta)
        if aliases[iv]:
            meta["aliases"] = aliases[iv]
        out.append(Seed(name, s.source, s.assembly, s.contig, s.start, s.end, meta))
    return out
