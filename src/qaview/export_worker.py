"""Child process for a figure export: `python -m qaview.export_worker <state.json> <out>`."""

import sys

from qaview.state import ViewState


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 2:
        print("usage: python -m qaview.export_worker <state.json> <output>", file=sys.stderr)
        return 2
    with open(argv[0]) as fh:
        state = ViewState.from_json(fh.read())
    print(f"[export] {state.describe()}", flush=True)
    state.render(argv[1])
    print(f"[export] wrote {argv[1]}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
