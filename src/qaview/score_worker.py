"""Child process for a scoring run:
`python -m qaview.score_worker <params.json> <recruitment.pkl> <out.pkl>`.

Kept separate from `score.py` so the parent never imports the pair-HMM engine or
compiles numba kernels on the server thread.
"""

import json
import sys

from qaview.score import run_scoring


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 3:
        print("usage: python -m qaview.score_worker <params.json> <recruitment.pkl> <out.pkl>",
              file=sys.stderr)
        return 2
    with open(argv[0]) as fh:
        params = json.load(fh)
    run_scoring(params, argv[1], argv[2])
    return 0


if __name__ == "__main__":
    sys.exit(main())
