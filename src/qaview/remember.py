"""Remembered alignments: snapshots of the right sidebar, kept in the output folder.

"Remember" saves what the panel shows -- the two stretches, the settings they
were made with, the alignment and, if they were shown, the reads -- as
`<output>/remembered/<id>.json`, plus `<id>.txt`, the alignment in blocks to
read without the app. Restoring loads the snapshot as saved, so it does not
change when thresholds or code do.

A snapshot is filed under the seeds of the locus it was made on but located by
genome coordinates, so its rectangle is drawn on any locus view whose segments
contain both stretches -- after re-recruitment, and for every seed that
recruits the same locus.
"""

from __future__ import annotations

import datetime
import json
import os
import re

REQUIRED = ("top", "mid", "bottom", "strand", "top_name", "top_contig", "top_start",
            "top_end", "bottom_name", "bottom_contig", "bottom_start", "bottom_end")
MAX_NOTE = 500
BLOCK = 60


def folder(output_dir: str) -> str:
    return os.path.join(output_dir, "remembered")


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)[:80]


def _meta(entry: dict) -> dict:
    """The entry without its alignment: what lists and rectangles need."""
    return {k: v for k, v in entry.items() if k != "alignment"}


def save(output_dir: str, alignment: dict, *, note: str, seeds: list, locus: str,
         settings: dict) -> dict:
    missing = [k for k in REQUIRED if k not in alignment]
    if missing:
        raise ValueError(f"not an alignment: missing {', '.join(missing)}")
    note = (note or "").strip()[:MAX_NOTE]
    now = datetime.datetime.now()
    ident = f"{_safe(seeds[0] if seeds else locus)}_{now:%Y%m%d-%H%M%S}"
    os.makedirs(folder(output_dir), exist_ok=True)
    while os.path.exists(os.path.join(folder(output_dir), f"{ident}.json")):
        ident += "_"
    a = alignment
    entry = {
        "id": ident, "created": now.isoformat(timespec="seconds"), "note": note,
        "seeds": list(seeds), "locus": locus, "settings": settings,
        "top": {"assembly": a["top_name"], "contig": a["top_contig"],
                "start": int(a["top_start"]), "end": int(a["top_end"])},
        "bottom": {"assembly": a["bottom_name"], "contig": a["bottom_contig"],
                   "start": int(a["bottom_start"]), "end": int(a["bottom_end"]),
                   "strand": a["strand"]},
        "reads": "reads" in a,
        "alignment": a,
    }
    base = os.path.join(folder(output_dir), ident)
    for ext, text in ((".json", json.dumps(entry)), (".txt", text_alignment(entry))):
        tmp = f"{base}{ext}.tmp"
        with open(tmp, "w") as fh:
            fh.write(text)
        os.replace(tmp, base + ext)
    return _meta(entry)


def load(output_dir: str, ident: str) -> dict:
    path = os.path.join(folder(output_dir), f"{_safe(ident)}.json")
    with open(path) as fh:
        return json.load(fh)


def forget(output_dir: str, ident: str) -> None:
    base = os.path.join(folder(output_dir), _safe(ident))
    if not os.path.exists(base + ".json"):
        raise FileNotFoundError(ident)
    for ext in (".json", ".txt"):
        if os.path.exists(base + ext):
            os.remove(base + ext)


def all_meta(output_dir: str) -> list[dict]:
    d = folder(output_dir)
    if not os.path.isdir(d):
        return []
    out = []
    for name in sorted(os.listdir(d)):
        if name.endswith(".json"):
            try:
                with open(os.path.join(d, name)) as fh:
                    out.append(_meta(json.load(fh)))
            except (OSError, ValueError):
                continue                                  # a torn file is skipped, not fatal
    return out


def place(metas: list[dict], segs_info1: list, segs_info2: list, asm1: str, asm2: str,
          flips1=None, flips2=None) -> list:
    """Plot rectangles for the snapshots whose stretches lie in this locus's layout.

    asm2 is the x axis and asm1 the y axis; `segs_info` rows are (contig, start,
    end, offset), `flips` flag segments drawn backwards. A stretch overlapping a
    segment is clipped to it.
    """
    def span(segs_info, flips, iv):
        for k, (contig, s, e, off) in enumerate(segs_info):
            if contig == iv["contig"] and iv["start"] < e and iv["end"] > s:
                lo, hi = max(iv["start"], s), min(iv["end"], e)
                if flips and flips[k]:
                    return off + e - hi, off + e - lo
                return off + lo - s, off + hi - s
        return None

    out = []
    for m in metas:
        if m["top"]["assembly"] != asm1 or m["bottom"]["assembly"] != asm2:
            continue
        y, x = span(segs_info1, flips1, m["top"]), span(segs_info2, flips2, m["bottom"])
        if x and y:
            out.append({"id": m["id"], "note": m["note"], "created": m["created"],
                        "x0": x[0], "x1": x[1], "y0": y[0], "y1": y[1]})
    return out


def text_alignment(entry: dict) -> str:
    """The alignment in blocks of BLOCK columns, 1-based coordinates at each line start."""
    a = entry["alignment"]
    t, b = entry["top"], entry["bottom"]
    rev = b["strand"] == "-"
    lines = [
        f"# qaview remembered alignment {entry['id']}",
        f"# created  {entry['created']}",
        f"# note     {entry['note']}" if entry["note"] else "# note     -",
        f"# seeds    {', '.join(entry['seeds'])}",
        f"# top      {t['assembly']} {t['contig']}:{t['start'] + 1:,}-{t['end']:,}",
        f"# bottom   {b['assembly']} {b['contig']}:{b['start'] + 1:,}-{b['end']:,}"
        + ("  (reverse complement)" if rev else ""),
        "# settings " + " ".join(f"{k}={v}" for k, v in entry["settings"].items()),
    ]
    if "compared" in a:
        lines.append(
            f"# compared {a['compared']:,} columns, identity {100 * a['identity']:.2f}%, "
            f"{a['mismatches']} mismatches, {a['ins']} bp in {a['ins_opens']} insertions, "
            f"{a['del']} bp in {a['del_opens']} deletions; outside the other stretch: "
            f"{a['outside_top']} bp top, {a['outside_bottom']} bp bottom")
    lines.append("# columns: '|' match, '.' mismatch, '-' gap, blank = outside the other stretch")
    lines.append("")
    w = max(len(t["assembly"]), len(b["assembly"]))
    nt, nb = 0, 0                                         # bases shown so far, per row
    for k in range(0, len(a["mid"]), BLOCK):
        top, mid, bot = a["top"][k:k + BLOCK], a["mid"][k:k + BLOCK], a["bottom"][k:k + BLOCK]
        pt = t["start"] + nt + 1
        pb = (b["end"] - nb) if rev else (b["start"] + nb + 1)
        # the position of the line's first base; none when the row has no base in it
        at_t = f"{pt:>13,}" if top.strip(" -") else " " * 13
        at_b = f"{pb:>13,}" if bot.strip(" -") else " " * 13
        lines.append(f"{t['assembly']:<{w}}  {at_t}  {top}")
        lines.append(f"{'':<{w}}  {'':>13}  {mid}")
        lines.append(f"{b['assembly']:<{w}}  {at_b}  {bot}")
        lines.append("")
        nt += sum(c not in " -" for c in top)
        nb += sum(c not in " -" for c in bot)
    return "\n".join(lines) + "\n"
