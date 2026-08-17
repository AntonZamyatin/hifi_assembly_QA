"""Child process for a recruitment run: `python -m qaview.recruit_worker <params.json> <out.pkl>`.

Kept separate from `recruit.py` so the parent never imports pysam-heavy recruiter state.
"""

import json
import sys

from qaview.recruit import run_recruitment


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 2:
        print("usage: python -m qaview.recruit_worker <params.json> <out.pkl>", file=sys.stderr)
        return 2
    with open(argv[0]) as fh:
        params = json.load(fh)
    run_recruitment(params, argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
