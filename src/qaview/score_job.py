"""Score one recruited locus into the project's score cache:
`python -m qaview.score_job <project.toml> <recruitment.pkl>`.

A subprocess so the server never imports the pair-HMM engine or forks from a
threaded process; this one forks the scoring workers itself. It goes through
`phmm_likelihood_report(cache_dir=...)`, the notebooks' entry point, so the key
and the lock are theirs and a locus is scored once whoever gets there first.
Prints `scored <i>/<n>` lines for the server's progress bar.
"""

import pickle
import sys


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 2:
        print("usage: python -m qaview.score_job <project.toml> <recruitment.pkl>",
              file=sys.stderr)
        return 2
    from qatool.plotters.likelihood import phmm_likelihood_report
    from qaview.project import load_project

    p = load_project(argv[0], check_files=False)
    with open(argv[1], "rb") as fh:
        rec = pickle.load(fh)
    step = max(1, len(rec["reads"]) // 100)
    last = [0]

    def progress(i, total):
        # Parallel scoring reports whole chunks, so print on progress, not on i % step.
        if i - last[0] >= step or i == total:
            last[0] = i
            print(f"scored {i}/{total}", flush=True)

    phmm_likelihood_report(
        rec["segments"], rec["reads"], p.asm1.name, p.asm2.name, db_path=p.db,
        bams={p.asm1.name: p.asm1.bam, p.asm2.name: p.asm2.bam},
        fastas={p.asm1.name: p.asm1.fasta, p.asm2.name: p.asm2.fasta},
        genome_fasta=p.genome_fasta, library_assembly=p.library_assembly,
        cache_dir=p.loci_dir, progress=progress)
    return 0


if __name__ == "__main__":
    sys.exit(main())
