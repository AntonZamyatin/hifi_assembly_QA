"""Regression test: the pair-HMM score on the HG002 seed locus.

Runs end to end on real data -- recruitment cache, aln database, BAMs, bgzipped
FASTAs -- through qatool's entry point, and skips when any of it is absent. The
locus is the one where hifiasm_hic and curated differ by a single event: an
(AT)n microsatellite three units shorter in hifiasm_hic, against the poly-T tail
of an Alu at h2tg000020l:61,904,626. The expectations are about that event, not
about a remembered number: which reads carry the verdict, and which way it points.

Run: venv/bin/python tests/test_score_hg002.py
"""

import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, ".qaview_cache", "recruit_3ca8364d592cb787.pkl")
DB = os.path.join(ROOT, "examples/dbs/hg002.aln.sqlite")
ASM = ("hifiasm_hic", "curated")
BAMS = {a: os.path.join(ROOT, f"examples/alignments/hg002.{a}.sorted.bam") for a in ASM}
FAS = {a: os.path.join(ROOT, f"examples/bgz/hg002.{a}.fasta.bgz") for a in ASM}

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def main():
    needed = [CACHE, DB, *BAMS.values(), *FAS.values()]
    missing = [p for p in needed if not os.path.exists(p)]
    if missing:
        print(f"SKIP: {len(missing)} data file(s) absent, e.g. {missing[0]}")
        return 0

    from qatool.plotters import phmm_likelihood_report

    d = pickle.load(open(CACHE, "rb"))
    r = phmm_likelihood_report(d["segments"], d["reads"], *ASM, db_path=DB,
                               bams=BAMS, fastas=FAS, genome_fasta=FAS["curated"],
                               library_assembly="curated")
    print()
    strong = [x for x in r.reads if x.ratio < -10]
    check(len(r.reads) == 106 and not r.dropped and not r.unscorable,
          f"expected all 106 recruited reads scored, got {len(r.reads)} "
          f"(+{len(r.dropped)} dropped, {len(r.unscorable)} unscorable)")
    check(r.delta_s < -200, f"expected a clear verdict for curated, got {r.delta_s:+.1f}")
    check(len(strong) >= 12, f"expected a dozen or more reads strongly for curated, got {len(strong)}")
    check(max(x.ratio for x in r.reads) < 1.0,
          "no read should favour hifiasm_hic: it carries the shorter allele")
    check(abs(r.exposure.value) < 0.1,
          f"segment lengths differ by 7 bp, so exposure should be ~0, got {r.exposure.value:+.3f}")
    check(r.worst_band_edge < 1e-6, f"band edge {r.worst_band_edge:.1e}: a band was too narrow")
    check(abs(r.at(50) - r.at(1000)) < 1e-6,
          "no read reaches the clamp here, so the verdict must not depend on C")

    # The speed-ups must not change the answer. Workers split the reads, so a
    # parallel run is the serial one reassembled; pruning skips placements only
    # after the boundary rule has seen all of them, so the read set is unchanged.
    from qamodel import fasta_length, library_stats, score_locus

    lib = library_stats(DB, "curated", fasta_length(FAS["curated"]))
    kw = dict(db_path=DB, bams=BAMS, fastas=FAS, library=lib)
    serial = score_locus(d["reads"], d["segments"], *ASM, workers=1, **kw)
    check(serial.reads == r.reads, "a parallel run must equal the serial one read for read")
    full = score_locus(d["reads"], d["segments"], *ASM, prune=None, workers=4, **kw)
    check([x.name for x in full.reads] == [x.name for x in r.reads]
          and full.dropped == r.dropped, "pruning must not change which reads are scored")
    drift = max(abs(a.ratio - b.ratio) for a, b in zip(full.reads, r.reads))
    check(drift < 1e-6, f"pruning moved a read's log-ratio by {drift:.2e}")

    print(f"all checks passed" if not _failures else f"{len(_failures)} FAILURE(S)")
    return 1 if _failures else 0


if __name__ == "__main__":
    sys.exit(main())
