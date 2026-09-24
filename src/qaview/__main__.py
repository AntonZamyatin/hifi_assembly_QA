"""Serve qaview: `python -m qaview [project.toml] [--port 5006]`.

Without a project file the page opens on a form to choose one (by its path on
this machine, with recently opened projects one click away). Opening a project
starts recruiting every seed and scoring every locus in the background
(docs/qaview.md §13.3); the page fills in as results arrive. Over SSH, forward
the port:

    ssh -L 5006:localhost:5006 <host>
    venv/bin/python -m qaview examples/projects/hg002_hic.toml

The earlier Panel app is still available as `panel serve src/qaview/app.py`.
"""

import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", nargs="?", help="project file (TOML); omit to choose one in the page")
    ap.add_argument("--port", type=int, default=5006)
    ap.add_argument("--host", default="127.0.0.1",
                    help="interface to listen on; keep the default and forward the port")
    args = ap.parse_args(argv)

    from qaview.project import ProjectError, load_project

    project = None
    if args.project:
        try:
            project = load_project(args.project)
        except ProjectError as exc:
            print(f"cannot open {args.project}:\n{exc}", file=sys.stderr)
            return 2
        print(f"{project.name}: {project.asm1.name} vs {project.asm2.name}, "
              f"{len(project.seeds)} seeds -> {project.output_dir}")
    else:
        print(f"no project given: open one in the page at http://{args.host}:{args.port}/")

    import uvicorn

    from qaview.server import create_app

    uvicorn.run(create_app(project), host=args.host, port=args.port, log_level="warning")
    return 0


if __name__ == "__main__":
    sys.exit(main())
