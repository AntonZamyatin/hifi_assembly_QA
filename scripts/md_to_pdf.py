#!/usr/bin/env python
"""Render a Markdown subset to PDF using matplotlib.

Exists because the mambaforge TeX installation on this machine has no format
files (`mktexfmt` is broken, `mktexlsr.pl` missing) and groff has no pdf device,
while repairing either would mean writing into the shared mambaforge base
environment.  matplotlib is already a project dependency, so this route needs no
install at all.

Supported: ATX headings, paragraphs, `-` bullets, fenced code blocks, `**bold**`,
`*italic*` and `` `code` `` inline spans, `---` rules, and `![caption](file.png)`
images.  Tables are not supported; write them as fenced code blocks so the
columns stay aligned.

Usage: venv/bin/python scripts/md_to_pdf.py in.md out.pdf [--title "..."]
"""

import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Arc, Circle, FancyArrowPatch

PAGE_W, PAGE_H = 8.27, 11.69          # A4 inches
MARGIN_L, MARGIN_R = 0.95, 0.95
MARGIN_T, MARGIN_B = 0.95, 0.85
COL_W = PAGE_W - MARGIN_L - MARGIN_R

SERIF = "DejaVu Serif"
MONO = "DejaVu Sans Mono"

BODY = 9.5
CODE = 8.0
CAPTION = 8.5
H_SIZE = {1: 17.0, 2: 13.0, 3: 10.8}
LEAD = 1.45                            # line height multiplier
PT = 1.0 / 72.0

# Bold alternative must come first so ** wins over *.
_INLINE = re.compile(r"(\*\*.+?\*\*|\*[^*\n]+?\*|`[^`]+?`)")


class Layout:
    """Greedy line-breaker that writes onto successive matplotlib pages."""

    def __init__(self, pdf, title=None):
        self.pdf = pdf
        self.title = title
        self.page = None
        self.page_no = 0
        self._measure_fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=100)
        self._renderer = self._measure_fig.canvas.get_renderer()
        self._new_page()

    # -- page handling ----------------------------------------------------
    def _new_page(self):
        if self.page is not None:
            self._footer()
            self.pdf.savefig(self.page)
            plt.close(self.page)
        self.page = plt.figure(figsize=(PAGE_W, PAGE_H))
        self.page_no += 1
        self.y = PAGE_H - MARGIN_T

    def _footer(self):
        self.page.text(PAGE_W / 2, MARGIN_B * 0.55, str(self.page_no),
                       transform=self._fig_tr(), ha="center", va="center",
                       fontname=SERIF, fontsize=8, color="0.45")

    def _fig_tr(self):
        # Inches -> figure fraction.
        return self.page.transFigure + _Scale(self.page)

    def close(self):
        self._footer()
        self.pdf.savefig(self.page)
        plt.close(self.page)
        plt.close(self._measure_fig)

    def space(self, inches):
        self.y -= inches

    def need(self, inches):
        if self.y - inches < MARGIN_B:
            self._new_page()

    # -- measurement ------------------------------------------------------
    def width_of(self, text, font, size, weight="normal", style="normal"):
        if not text:
            return 0.0
        t = self._measure_fig.text(0, 0, text, fontname=font, fontsize=size,
                                   fontweight=weight, fontstyle=style)
        w = t.get_window_extent(self._renderer).width / self._measure_fig.dpi
        t.remove()
        return w

    # -- drawing ----------------------------------------------------------
    def draw(self, x_in, text, font, size, weight="normal", style="normal",
             color="black"):
        self.page.text(x_in / PAGE_W, self.y / PAGE_H, text,
                       fontname=font, fontsize=size, fontweight=weight,
                       fontstyle=style, color=color, ha="left", va="top")

    def runs(self, text):
        """Split inline markup into (text, font, weight, style) runs."""
        out = []
        for piece in _INLINE.split(text):
            if not piece:
                continue
            if piece.startswith("**") and piece.endswith("**"):
                out.append((piece[2:-2], SERIF, "bold", "normal"))
            elif piece.startswith("`") and piece.endswith("`"):
                out.append((piece[1:-1], MONO, "normal", "normal"))
            elif len(piece) > 2 and piece.startswith("*") and piece.endswith("*"):
                out.append((piece[1:-1], SERIF, "normal", "italic"))
            else:
                out.append((piece, SERIF, "normal", "normal"))
        return out

    def tokens(self, text):
        """(word, font, weight, space_before).

        The flag matters: splitting runs naively puts a space before any
        punctuation that follows an inline span, so `code`. renders as "code ."
        """
        out, pending = [], False
        for run_text, font, weight, style in self.runs(text):
            parts = [w for w in run_text.split(" ") if w]
            if run_text.startswith(" "):
                pending = True
            for j, w in enumerate(parts):
                out.append((w, font, weight, style, pending or j > 0))
                pending = False
            if run_text.endswith(" ") or not parts:
                pending = True
        if out:
            w, f, wt, st, _ = out[0]
            out[0] = (w, f, wt, st, False)
        return out

    def paragraph(self, text, size=BODY, indent=0.0, hang=0.0, bullet=None):
        line_h = size * LEAD * PT
        avail = COL_W - indent
        space_w = self.width_of(" ", SERIF, size)

        line, line_w, first = [], 0.0, True
        for w, font, weight, style, sp in self.tokens(text):
            ww = self.width_of(w, font, size, weight, style)
            gap = space_w if (sp and line) else 0.0
            limit = avail - (0.0 if first else hang)
            if line and line_w + gap + ww > limit:
                self._flush(line, size, indent, 0.0 if first else hang,
                            bullet if first else None, space_w)
                self.y -= line_h
                self.need(line_h)
                line, line_w, first = [(w, font, weight, style, False)], ww, False
            else:
                line.append((w, font, weight, style, sp and bool(line)))
                line_w += gap + ww
        if line:
            self._flush(line, size, indent, 0.0 if first else hang,
                        bullet if first else None, space_w)
            self.y -= line_h
        self.need(line_h)

    def _flush(self, line, size, indent, hang, bullet, space_w):
        x = MARGIN_L + indent + hang
        if bullet is not None:
            self.draw(MARGIN_L + indent - 0.16, bullet, SERIF, size)
        for w, font, weight, style, sp in line:
            if sp:
                x += space_w
            self.draw(x, w, font, size, weight, style)
            x += self.width_of(w, font, size, weight, style)

    def heading(self, level, text):
        size = H_SIZE[level]
        line_h = size * LEAD * PT
        self.space(0.16 if level > 1 else 0.10)
        self.need(line_h * 3)
        if level == 1 and self.page_no > 1 and self.y < PAGE_H - MARGIN_T - 0.01:
            self._new_page()
        self.draw(MARGIN_L, text, SERIF, size, "bold")
        self.y -= line_h
        if level == 1:
            self.y -= 0.04
            self.page.add_artist(plt.Line2D(
                [MARGIN_L / PAGE_W, (PAGE_W - MARGIN_R) / PAGE_W],
                [self.y / PAGE_H, self.y / PAGE_H],
                transform=self.page.transFigure, color="0.75", linewidth=0.8))
            self.y -= 0.10
        self.space(0.04)

    def code_block(self, lines):
        line_h = CODE * LEAD * PT
        block_h = line_h * len(lines) + 0.16
        self.space(0.06)
        if block_h < (PAGE_H - MARGIN_T - MARGIN_B):
            self.need(block_h)
        top = self.y
        self.page.add_artist(plt.Rectangle(
            ((MARGIN_L - 0.10) / PAGE_W, (top - block_h + 0.06) / PAGE_H),
            (COL_W + 0.20) / PAGE_W, (block_h - 0.04) / PAGE_H,
            transform=self.page.transFigure, facecolor="0.955",
            edgecolor="0.88", linewidth=0.6, zorder=0))
        self.y -= 0.08
        for ln in lines:
            self.need(line_h)
            self.draw(MARGIN_L, ln.rstrip(), MONO, CODE)
            self.y -= line_h
        self.space(0.12)

    def image(self, path, caption=None):
        """Place a raster figure at column width, with an optional caption."""
        img = plt.imread(path)
        h_in = COL_W * img.shape[0] / img.shape[1]
        cap_h = 0.20 if caption else 0.0
        self.space(0.10)
        self.need(h_in + cap_h + 0.10)
        ax = self.page.add_axes([MARGIN_L / PAGE_W, (self.y - h_in) / PAGE_H,
                                 COL_W / PAGE_W, h_in / PAGE_H])
        ax.imshow(img)
        ax.axis("off")
        self.y -= h_in + 0.04
        if caption:
            self.draw(MARGIN_L, caption, SERIF, CAPTION, style="italic", color="0.35")
            self.y -= cap_h
        self.space(0.08)

    def rule(self):
        self.space(0.08)
        self.need(0.1)
        self.page.add_artist(plt.Line2D(
            [MARGIN_L / PAGE_W, (PAGE_W - MARGIN_R) / PAGE_W],
            [self.y / PAGE_H, self.y / PAGE_H],
            transform=self.page.transFigure, color="0.8", linewidth=0.7))
        self.space(0.12)


class _Scale:
    """Placeholder so _fig_tr stays valid; footer uses figure fractions."""

    def __init__(self, fig):
        self.fig = fig

    def __radd__(self, other):
        return other


def state_diagram(pdf):
    """Draw the transition graph as a real figure rather than ASCII.

    The two forbidden edges are named in the caption rather than drawn dashed:
    routing them through the graph crosses every other edge and makes the
    picture unreadable, which defeats the point of having one.
    """
    fig = plt.figure(figsize=(PAGE_W, 5.2))
    ax = fig.add_axes([0.03, 0.02, 0.94, 0.84])
    ax.set_xlim(0, 10)
    ax.set_ylim(-2.1, 5.6)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.text(0.5, 0.955, "The pair HMM transition graph",
             ha="center", fontname=SERIF, fontsize=12, fontweight="bold")
    fig.text(0.5, 0.915,
             "Two edges are deliberately absent: B→D and D→E. "
             "Both remove redundancy, not real possibilities.",
             ha="center", fontname=SERIF, fontsize=9, color="0.3")

    pos = {"B": (0.9, 2.0), "M": (3.4, 2.0), "I": (5.9, 3.9),
           "D": (5.9, 0.1), "E": (9.0, 2.0)}
    R = 0.34
    for name, (x, y) in pos.items():
        ax.add_patch(Circle((x, y), R, facecolor="white",
                            edgecolor="0.2", linewidth=1.3, zorder=3))
        ax.text(x, y, f"${name}$", ha="center", va="center", zorder=4,
                fontname=SERIF, fontsize=13)

    def edge(a, b, label, rad=0.0, nudge=(0.0, 0.0)):
        (x1, y1), (x2, y2) = pos[a], pos[b]
        ax.add_patch(FancyArrowPatch(
            (x1, y1), (x2, y2), connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>", mutation_scale=12,
            shrinkA=R * 72 * 0.42, shrinkB=R * 72 * 0.42,
            linewidth=1.1, color="0.3", zorder=2))
        # matplotlib's Arc3 control point is mid + rad*(dy, -dx); the quadratic
        # Bezier midpoint is therefore mid + (rad/2)*(dy, -dx). Computing it
        # beats hand-placing every label and then chasing them when a node moves.
        dx, dy = x2 - x1, y2 - y1
        lx = (x1 + x2) / 2 + rad / 2 * dy + nudge[0]
        ly = (y1 + y2) / 2 - rad / 2 * dx + nudge[1]
        ax.text(lx, ly, label, ha="center", va="center", fontname=SERIF,
                fontsize=9, color="0.15", zorder=5,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="none", alpha=0.95))

    edge("B", "M", r"$1-\delta_I$", 0.00, (0.0, 0.18))
    edge("B", "I", r"$\delta_I$", -0.28)
    edge("M", "I", r"$\delta_I$", 0.26)
    edge("I", "M", r"$1-\gamma_I-\tau$", 0.26)
    edge("M", "D", r"$\delta_D$", 0.26)
    edge("D", "M", r"$1-\gamma_D$", 0.26)
    edge("M", "E", r"$\tau$", 0.00, (0.0, 0.18))
    edge("I", "E", r"$\tau$", -0.24)

    def self_loop(state, above, label, side=None):
        x, y = pos[state]
        cy = y + (0.68 if above else -0.68)
        ax.add_patch(Arc((x, cy), 0.68, 0.64,
                         theta1=(-15 if above else 165),
                         theta2=(195 if above else 375),
                         edgecolor="0.3", linewidth=1.1, zorder=2))
        ax.plot([x + (0.33 if above else -0.33)],
                [cy + (-0.08 if above else 0.08)],
                marker=(3, 0, 20 if above else 200),
                markersize=6, color="0.3", zorder=3)
        if side == "left":
            ax.text(x - 0.60, cy, label, ha="right", va="center",
                    fontname=SERIF, fontsize=9, color="0.15")
        else:
            ax.text(x, y + (1.34 if above else -1.34), label, ha="center",
                    va="center", fontname=SERIF, fontsize=9, color="0.15")

    self_loop("M", True, r"$1-\delta_I-\delta_D-\tau$", side="left")
    self_loop("I", True, r"$\gamma_I$")
    self_loop("D", False, r"$\gamma_D$")

    notes = ("M consumes one read base and one template base    "
             "I consumes a read base only    "
             "D consumes a template base only")
    ax.text(5.0, -1.85, notes, ha="center", va="center", fontname=SERIF,
            fontsize=8.6, color="0.4")

    pdf.savefig(fig)
    plt.close(fig)


def render(md_path, pdf_path, title=None, diagram_after=None):
    with open(md_path) as fh:
        lines = fh.read().split("\n")

    with PdfPages(pdf_path) as pdf:
        lay = Layout(pdf, title)
        i, para = 0, []

        def flush_para():
            if para:
                lay.paragraph(" ".join(para))
                lay.space(0.07)
                para.clear()

        while i < len(lines):
            ln = lines[i]
            stripped = ln.strip()

            if stripped.startswith("```"):
                flush_para()
                i += 1
                block = []
                while i < len(lines) and not lines[i].strip().startswith("```"):
                    block.append(lines[i])
                    i += 1
                lay.code_block(block)
                i += 1
                continue

            m = re.match(r"^!\[(.*?)\]\((.+?)\)$", stripped)
            if m:
                flush_para()
                lay.image(m.group(2), m.group(1) or None)
                i += 1
                continue

            m = re.match(r"^(#{1,3})\s+(.*)$", stripped)
            if m:
                flush_para()
                lay.heading(len(m.group(1)), m.group(2))
                if diagram_after and m.group(2).strip() == diagram_after:
                    lay.space(0.02)
                i += 1
                continue

            if stripped in ("---", "***", "___"):
                flush_para()
                lay.rule()
                i += 1
                continue

            if stripped.startswith("- "):
                flush_para()
                # A bullet runs until a blank line, the next bullet, or any
                # other block; wrapped source lines belong to the same bullet.
                item = [stripped[2:]]
                i += 1
                while i < len(lines):
                    nxt = lines[i].strip()
                    if (not nxt or nxt.startswith("- ") or nxt.startswith("#")
                            or nxt.startswith("```") or nxt.startswith("![")
                            or nxt in ("---", "***", "___")):
                        break
                    item.append(nxt)
                    i += 1
                lay.paragraph(" ".join(item), indent=0.26, hang=0.0, bullet="•")
                lay.space(0.03)
                continue

            if not stripped:
                flush_para()
                i += 1
                continue

            para.append(stripped)
            i += 1

        flush_para()
        lay.close()

        if diagram_after:
            state_diagram(pdf)

        info = pdf.infodict()
        info["Title"] = title or md_path
        info["Author"] = "Anton Zamyatin"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--title", default=None)
    ap.add_argument("--diagram-after", default=None,
                    help="append the pair-HMM state diagram as a final page")
    a = ap.parse_args()
    render(a.src, a.dst, a.title, a.diagram_after)
    print(f"wrote {a.dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
