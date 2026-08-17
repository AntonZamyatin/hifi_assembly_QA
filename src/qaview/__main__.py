"""Serve the qaview app: `python -m qaview [--port 5006] [--show]`.

Equivalent to `panel serve src/qaview/app.py`, but importable so the entry point does
not depend on the source layout.  Over SSH, forward the port rather than using --show:

    ssh -L 5006:localhost:5006 <host>
    python -m qaview --port 5006
"""

import argparse
import sys

import panel as pn

from qaview.app import build_app


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=5006)
    ap.add_argument("--show", action="store_true", help="open a browser (local sessions only)")
    ap.add_argument("--allow-websocket-origin", action="append", default=[],
                    metavar="HOST:PORT",
                    help="extra allowed origin; repeatable.  localhost and 127.0.0.1 "
                         "on --port are always allowed")
    ap.add_argument("--cache-dir", default=".qaview_cache",
                    help="where recruitment results are cached")
    args = ap.parse_args(argv)

    # Bokeh matches the browser's Origin header literally, so localhost and 127.0.0.1
    # are distinct even though they resolve to the same host.  Allow both, or the page
    # loads and then sits dead with a 403 on /ws.
    origins = [f"localhost:{args.port}", f"127.0.0.1:{args.port}"]
    origins += [o for o in args.allow_websocket_origin if o not in origins]

    pn.serve(lambda: build_app(cache_dir=args.cache_dir), port=args.port,
             show=args.show, title="qaview", allow_websocket_origin=origins)
    return 0


if __name__ == "__main__":
    sys.exit(main())
