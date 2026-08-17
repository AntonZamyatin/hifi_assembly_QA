#!/usr/bin/env python
"""Render an interactive qaview report to a standalone HTML file.

Recruits from a seed locus (or skips recruitment with --no-recruit) and writes a
self-contained page — bokeh JS is inlined, so it opens over scp/sshfs with no server.

    scripts/qaview_render.py --base examples --label rVipUrs1 \
        --contig SUPER_Z --start 5404885 --end 5409029 -o /tmp/rVipUrs1.html

Phase 1: static panels with linked zoom and a shared crosshair.  No widgets yet, so
every threshold is fixed at the values passed here.
"""

import argparse
import os
import sys

from qatool.SegmentStorage import Segment, SegmentCollection
from qaview.glyphs import report_layout, save_html
from qaview.region import RegionData


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="examples", help="directory holding dbs/ alignments/ bgz/")
    ap.add_argument("--label", required=True, help="dataset label, e.g. rVipUrs1")
    ap.add_argument("--asm1", default="ref")
    ap.add_argument("--asm2", default="lja")
    ap.add_argument("--contig", required=True, help="seed contig, in asm1")
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("-o", "--out", default="qaview.html")

    ap.add_argument("--no-recruit", action="store_true",
                    help="skip recruitment; needs --asm2-contig/--asm2-start/--asm2-end")
    ap.add_argument("--asm2-contig")
    ap.add_argument("--asm2-start", type=int)
    ap.add_argument("--asm2-end", type=int)

    ap.add_argument("--pid-threshold", type=float, default=0.998, help="for plotting")
    ap.add_argument("--recruit-pid", type=float, default=0.99)
    ap.add_argument("--plen-threshold", type=float, default=0.999)
    ap.add_argument("--best-pid-delta", type=float, default=0.001)
    ap.add_argument("--pid-tie-eps", type=float, default=0.001)
    ap.add_argument("--binsize", type=int, default=None,
                help="default: auto (~1000 bins on the longer axis)")
    ap.add_argument("--merge-gap", type=int, default=50_000)
    ap.add_argument("--ext-window", type=int, default=1000)
    ap.add_argument("--ext-max", type=int, default=300_000)
    ap.add_argument("--ext-rounds", type=int, default=2)
    ap.add_argument("--max-rounds", type=int, default=100)
    ap.add_argument("--recruit-out", default=None,
                    help="recruiter scratch dir (default: alongside --out)")
    ap.add_argument("--color-by", choices=("class", "block"), default="class")
    args = ap.parse_args(argv)

    b = args.base
    paths = dict(
        db=f"{b}/dbs/{args.label}.aln.sqlite",
        bam1=f"{b}/alignments/{args.label}.{args.asm1}.sorted.bam",
        bam2=f"{b}/alignments/{args.label}.{args.asm2}.sorted.bam",
        fa1=f"{b}/bgz/{args.label}.{args.asm1}.fasta.bgz",
        fa2=f"{b}/bgz/{args.label}.{args.asm2}.fasta.bgz",
    )
    missing = [p for p in paths.values() if not os.path.exists(p)]
    if missing:
        ap.error("missing input files:\n  " + "\n  ".join(missing))

    if args.no_recruit:
        if not all((args.asm2_contig, args.asm2_start is not None, args.asm2_end is not None)):
            ap.error("--no-recruit needs --asm2-contig, --asm2-start and --asm2-end")
        segments = SegmentCollection([
            Segment(args.asm1, args.contig, args.start, args.end),
            Segment(args.asm2, args.asm2_contig, args.asm2_start, args.asm2_end),
        ])
        reads: set = set()
    else:
        from qatool.SRRecruiter import Recruiter
        out_dir = args.recruit_out or os.path.join(
            os.path.dirname(os.path.abspath(args.out)) or ".", f"{args.label}_recruit")
        os.makedirs(out_dir, exist_ok=True)
        print(f"[qaview] recruiting from {args.contig}:{args.start:,}-{args.end:,} "
              f"(scratch: {out_dir})", flush=True)
        srr = Recruiter(
            db_path=paths["db"], asm1=args.asm1, asm2=args.asm2,
            asm1_bam_path=paths["bam1"], asm2_bam_path=paths["bam2"],
            asm1_fasta_path=paths["fa1"], asm2_fasta_path=paths["fa2"],
            merge_gap=args.merge_gap, pid_threshold=args.recruit_pid,
            plen_threshold=args.plen_threshold, best_pid_delta=args.best_pid_delta,
            extension_window=args.ext_window, ext_max=args.ext_max,
            ext_rounds=args.ext_rounds, max_rounds=args.max_rounds,
            out_dir=out_dir, verbose=True,
        )
        segments, reads, _df = srr.run(
            SegmentCollection([Segment(args.asm1, args.contig, args.start, args.end)]))

    print("[qaview] building RegionData ...", flush=True)
    rd = RegionData.build(
        segments, reads, paths["db"], paths["bam1"], paths["bam2"],
        asm1_fasta=paths["fa1"], asm2_fasta=paths["fa2"],
        asm1=args.asm1, asm2=args.asm2,
        pid_threshold=args.pid_threshold, plen_threshold=args.plen_threshold,
        binsize=args.binsize, best_pid_delta=args.best_pid_delta,
        pid_tie_eps=args.pid_tie_eps)

    total_blocks = len(rd.lines) + rd.lines.n_dropped
    pct = (100.0 * rd.lines.n_dropped / total_blocks) if total_blocks else 0.0
    print(f"[qaview] {len(rd.read_names)} reads, {len(rd.lines)} lines, "
          f"{rd.lines.n_dropped} blocks dropped at segment edges ({pct:.1f}%)")
    if pct > 25:
        print("[qaview] WARNING: a large share of blocks fall outside single segments, "
              "so the dotplot understates real correspondence (docs/qaview.md §6).")

    save_html(report_layout(rd, color_by=args.color_by), args.out,
              title=f"qaview {args.label} {args.contig}:{args.start}-{args.end}")
    print(f"[qaview] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
