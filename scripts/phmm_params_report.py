#!/usr/bin/env python
"""Build the stage-G parameter report: how alignments were selected, what events
were counted, and what the counts imply.

    venv/bin/python scripts/phmm_params_report.py \
        --json docs/phmm_params_hg002.json \
        --bam examples/alignments/hg002.curated.sorted.bam \
        --out docs/phmm_params_report.pdf

The headline statistics come from the JSON `qamodel.estimate` writes. The
between-tile figure needs per-tile counts, which that JSON does not carry, so
the script re-scans a smaller sample (~50 s) and caches it next to the output.
"""

import argparse
import json
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from qamodel.estimate import Counts, estimate, scan_bam       # noqa: E402

# Categorical slots 1 and 2 of the reference palette, validated as a pair
# (CVD dE 24.7, normal-vision 33.6, both >= 3:1 on this surface).
INS, DEL = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d8d8d4"
FIG_DPI = 200


def style(ax):
    """Recessive axes: the data is the ink, the frame is not."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=7.5, length=3)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for lab in (ax.xaxis.label, ax.yaxis.label):
        lab.set_color(MUTED)
        lab.set_fontsize(8)
    ax.title.set_color(INK)
    ax.title.set_fontsize(9)


def fig_run_lengths(counts, params, path):
    """Observed gap run lengths against the geometric law the model assumes.

    A gap of length n has probability (1-gamma) gamma^(n-1) under the model, so
    on a log axis the prediction is a straight line. Where the observed bars sit
    above it, the model cannot produce gaps that the data contains.
    """
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.9), dpi=FIG_DPI)
    for ax, hist, gamma, color, name in (
            (axes[0], counts["i_hist"], params["gamma_i"], INS, "insertion"),
            (axes[1], counts["d_hist"], params["gamma_d"], DEL, "deletion")):
        h = np.asarray(hist[1:], dtype=float)
        n = np.arange(1, len(h) + 1)
        obs = h / h.sum()
        pred = (1 - gamma) * gamma ** (n - 1)
        ax.bar(n, obs, color=color, width=0.72, label="observed")
        ax.plot(n, pred, color=INK, linewidth=1.6, linestyle=(0, (4, 2)),
                label=f"geometric, $\\gamma$ = {gamma:.3f}")
        ax.set_yscale("log")
        ax.set_ylim(1e-7, 2.0)
        ax.set_xlim(0.3, 12.5)
        ax.set_xticks(range(1, 13))
        ax.set_title(f"{name} run lengths")
        ax.set_xlabel("run length (bp)")
        style(ax)
        ax.legend(frameon=False, fontsize=7.5, labelcolor=MUTED)
    axes[0].set_ylabel("share of runs")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_calibration(cal, path):
    """Reported quality against the substitution rate actually observed."""
    bases = np.asarray(cal["bases"], dtype=float)
    mism = np.asarray(cal["mismatches"], dtype=float)
    q = np.flatnonzero(bases > 10_000)
    obs = mism[q] / bases[q]
    reported = 10.0 ** (-q / 10.0)

    fig, ax = plt.subplots(figsize=(4.9, 3.4), dpi=FIG_DPI)
    lim = [3e-6, 1.2]
    ax.plot(lim, lim, color=MUTED, linewidth=1.0, linestyle=(0, (4, 2)), zorder=1)
    ax.set_aspect("equal")            # so the diagonal really is at 45 degrees
    ax.text(1.2e-4, 1.9e-4, "reported = observed", fontsize=7.5, color=MUTED,
            rotation=45, rotation_mode="anchor", ha="left", va="bottom")
    sizes = 14 + 46 * (bases[q] / bases[q].max()) ** 0.5
    ax.scatter(reported, obs, s=sizes, color=INS, zorder=3,
               edgecolor="white", linewidth=1.2)
    for qq, r, o in zip(q, reported, obs):
        ax.annotate(f"Q{qq}", (r, o), textcoords="offset points",
                    xytext=(9, -9), fontsize=7.5, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("error probability implied by the reported quality")
    ax.set_ylabel("observed substitution rate")
    ax.set_title("Reported qualities are pessimistic at every level")
    style(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_context(cal, path):
    """Gap rate inside homopolymers against everywhere else."""
    hp = cal["hp_gaps"] / max(cal["hp_bases"], 1)
    other = cal["other_gaps"] / max(cal["other_bases"], 1)
    fig, ax = plt.subplots(figsize=(4.4, 2.5), dpi=FIG_DPI)
    # One hue: both bars are the same quantity measured in two contexts, and
    # the other figures already use blue/orange for insertion/deletion.
    bars = ax.barh([1, 0], [hp, other], color=INS, height=0.55)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["inside a homopolymer\nof 4 bp or more", "everywhere else"],
                       fontsize=8, color=INK)
    for bar, v in zip(bars, (hp, other)):
        ax.text(v + 1.4e-4, bar.get_y() + bar.get_height() / 2, f"{v:.2e} / bp",
                va="center", fontsize=8, color=INK)
    ax.set_xlim(0, max(hp, other) * 1.35)
    ax.set_xlabel("gap events per aligned base")
    ax.set_title(f"Gaps are {hp / other:.1f}x more frequent in homopolymers")
    style(ax)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_between_tiles(per_tile, jk, path):
    """Per-tile estimates of delta_I, against the binomial error bar.

    The point of the figure: the spread of the histogram is what the jackknife
    measures, and it is far wider than the binomial standard error, which
    assumes every aligned base is an independent draw.
    """
    vals = np.asarray([t["delta_i"] for t in per_tile]) * 1e4
    centre = np.median(vals)
    fig, ax = plt.subplots(figsize=(5.4, 2.9), dpi=FIG_DPI)
    ax.hist(vals, bins=np.linspace(0, 12, 70), color=INS)
    ax.axvline(centre, color=INK, linewidth=1.4)
    ax.annotate(f"pooled estimate {jk['delta_i_estimate'] * 1e4:.2f}",
                (centre, ax.get_ylim()[1] * 0.92), xytext=(10, 0),
                textcoords="offset points", fontsize=8, color=INK)
    ax.set_xlabel(r"per-tile $\delta_I$  ($\times 10^{-4}$)")
    ax.set_ylabel("tiles")
    ax.set_title("One tile is not one draw: the between-tile spread")
    style(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def per_tile_counts(bam, cache, n_windows=4000, max_reads=100_000):
    if os.path.exists(cache):
        with open(cache) as fh:
            return json.load(fh)
    print(f"scanning {n_windows:,} tiles for the between-tile figure ...")
    _, per = scan_bam(bam, n_windows=n_windows, max_reads=max_reads,
                      max_per_window=25, per_window=True)
    out = []
    for c in per:
        m_total = c.m_m + c.m_i + c.m_d
        if m_total < 50_000:                  # a tile with too little sequence
            continue
        p = estimate(c).params
        out.append({"delta_i": p.delta_i, "delta_d": p.delta_d,
                    "gamma_i": p.gamma_i, "gamma_d": p.gamma_d, "m": m_total})
    with open(cache, "w") as fh:
        json.dump(out, fh)
    return out


def table(rows, widths=None):
    """Fixed-width table as code-block lines, which is what md_to_pdf renders."""
    widths = widths or [max(len(str(r[i])) for r in rows) for i in range(len(rows[0]))]
    out = []
    for r in rows:
        cells = [str(c).ljust(w) if i == 0 else str(c).rjust(w)
                 for i, (c, w) in enumerate(zip(r, widths))]
        out.append("  ".join(cells).rstrip())
    return out


def build_markdown(d, per_tile, figs, bam):
    c, p, jk, cal = d["counts"], d["params"], d["jackknife"], d["calibration"]
    f = d["filters"]
    seen = c["reads"] + c["skip_mapq"] + c["skip_split_read"] + c["skip_clip"] \
        + c["skip_identity"]
    gap_events = c["m_i"] + c["m_d"]
    cal_bases, cal_mism = sum(cal["bases"]), sum(cal["mismatches"])
    mm_rate = cal_mism / cal_bases
    hp = cal["hp_gaps"] / max(cal["hp_bases"], 1)
    other = cal["other_gaps"] / max(cal["other_bases"], 1)

    def pct(x, of):
        return f"{100 * x / of:.2f}%"

    sel = table([
        ["filter", "kept", "dropped", "share"],
        ["primary alignments in the sampled tiles", f"{seen:,}", "", ""],
        ["mapping quality >= " + str(f["min_mapq"]), "", f"{c['skip_mapq']:,}",
         pct(c["skip_mapq"], seen)],
        ["no SA tag (not a split read)", "", f"{c['skip_split_read']:,}",
         pct(c["skip_split_read"], seen)],
        ["clipped fraction <= 5%", "", f"{c['skip_clip']:,}", pct(c["skip_clip"], seen)],
        [f"divergence de <= {f['max_de']}", "", f"{c['skip_identity']:,}",
         pct(c["skip_identity"], seen)],
        ["counted", f"{c['reads']:,}", "", pct(c["reads"], seen)],
    ])

    ev = table([
        ["event", "count", "per aligned base"],
        ["aligned columns (M)", f"{c['m_bases']:,}", ""],
        ["  matches", f"{c['m_bases'] - c['mismatches']:,}", f"{1 - c['mismatches'] / c['m_bases']:.6f}"],
        ["  mismatches (NM)", f"{c['mismatches']:,}", f"{c['mismatches'] / c['m_bases']:.3e}"],
        ["  mismatches (base by base)", f"{cal_mism:,}", f"{mm_rate:.3e}"],
        ["insertion opens (M->I)", f"{c['m_i']:,}", f"{c['m_i'] / c['m_bases']:.3e}"],
        ["insertion extends (I->I)", f"{c['i_i']:,}", f"{c['i_i'] / c['m_bases']:.3e}"],
        ["insertion closes (I->M)", f"{c['i_m']:,}", ""],
        ["inserted bases", f"{c['i_bases']:,}", f"{c['i_bases'] / c['m_bases']:.3e}"],
        ["deletion opens (M->D)", f"{c['m_d']:,}", f"{c['m_d'] / c['m_bases']:.3e}"],
        ["deletion extends (D->D)", f"{c['d_d']:,}", f"{c['d_d'] / c['m_bases']:.3e}"],
        ["deletion closes (D->M)", f"{c['d_m']:,}", ""],
        ["deleted bases", f"{c['d_bases']:,}", f"{c['d_bases'] / c['m_bases']:.3e}"],
        ["gap events of either kind", f"{gap_events:,}", f"{gap_events / c['m_bases']:.3e}"],
        ["paths split at a run > %d bp" % f["max_gap"], f"{c['split']:,}", ""],
        ["I<->D adjacencies (no model edge)", f"{c['i_d'] + c['d_i']:,}", ""],
    ])

    est = table([
        ["parameter", "estimate", "jackknife", "binomial", "meaning"],
        ["delta_I", f"{p['delta_i']:.4e}", f"{jk['delta_i']:.1e}",
         f"{d['se']['delta_i']:.1e}", "insertion opens per M column"],
        ["delta_D", f"{p['delta_d']:.4e}", f"{jk['delta_d']:.1e}",
         f"{d['se']['delta_d']:.1e}", "deletion opens per M column"],
        ["gamma_I", f"{p['gamma_i']:.4f}", f"{jk['gamma_i']:.4f}",
         f"{d['se']['gamma_i']:.4f}", f"mean run {1 / (1 - p['gamma_i']):.2f} bp"],
        ["gamma_D", f"{p['gamma_d']:.4f}", f"{jk['gamma_d']:.4f}",
         f"{d['se']['gamma_d']:.4f}", f"mean run {1 / (1 - p['gamma_d']):.2f} bp"],
    ])

    qrows = [["reported", "bases", "implied eps", "observed", "as Phred"]]
    bases = np.asarray(cal["bases"], dtype=float)
    for q in np.flatnonzero(bases > 100_000):
        o = cal["mismatches"][q] / bases[q]
        qrows.append([f"Q{q}", f"{int(bases[q]):,}", f"{10 ** (-q / 10):.2e}",
                      f"{o:.2e}", f"Q{-10 * np.log10(o):.1f}" if o else "-"])

    return f"""# Pair-HMM parameters from HiFi alignments

Stage G of `docs/phmm.md`: the four gap parameters of the pair HMM, estimated
from reads aligned to the HG002 curated assembly. This report documents which
alignments were used, what was counted in them, and what the counts imply.

Source: `{os.path.basename(bam)}`, minimap2 2.30 `map-hifi`. Generated by
`qamodel/estimate.py`; provenance in `docs/phmm_params_hg002.json`.

## 1. Choosing where to look

The assembly is tiled into non-overlapping {f['windows'] and 100_000 // 1000}
kb tiles and {f['windows']:,} of them are drawn **without replacement** from a
seeded RNG, so a run reproduces from `--seed {f['seed']}` alone. That is
{f['windows'] * 100_000 / 1e9:.1f} Gb of the 6.0 Gb diploid assembly.

Sampling tiles rather than reading the head of the BAM is the point: gap rates
vary several-fold with sequence context, and context is not uniformly
distributed along an assembly. Reading the first N reads would give one region's
opinion at full apparent precision.

The draws are without replacement because overlapping windows would count a read
in the intersection twice, and would break the jackknife of section 5, which
treats each tile as an independent block.

## 2. Choosing which alignments to count

Within a tile, every alignment touching it is offered, and rejected in this
order. Secondary, supplementary and unmapped records are excluded first —
secondary records carry no sequence or qualities at all. An alignment whose
start lies before the tile was counted in an earlier tile and is skipped, so
each read is attributed to exactly one tile.

```
{chr(10).join(sel)}
```

At most {25} reads are taken from any one tile, so the budget spreads over many
tiles instead of exhausting the first few.

**Mapping quality is not a uniqueness filter here.** The target is a *diploid*
assembly: a read from a homozygous region maps equally well to both haplotypes
and scores MAPQ 0 however perfect its alignment. Requiring MAPQ 60 keeps under
10% of primaries — exactly the haplotype-divergent ones. Measured both ways the
gap estimates agree within about 3 standard errors, so the threshold is 0 and
the identity, clipping and split-read filters do the work.

Each remaining filter removes alignments whose *path* cannot be believed: a
split read is least reliable exactly at its ends, a heavily clipped read is
mostly somewhere else, and a read at 3% divergence is evidence about the
assembly rather than about sequencing error.

## 3. What counts as an event

The model's states are M, I and D, so an alignment is already a state path:
`100M1I100M` is 100 M columns, one I, then 100 more. `M`, `=` and `X` all map to
M — the state does not care whether the base matched, only that read and
template both advanced. A run of length n contributes n-1 internal transitions
plus one transition into it from the run before.

- **Clips, reference skips and padding break the path** rather than
  transitioning through it. A clip means the read left the template, which the
  model's free start and free end already handle.
- **Runs longer than {f['max_gap']} bp split the path.** A 138 bp deletion is a
  variant, an assembly error or a collapsed repeat, not a sequencing indel, and
  counting it drives gamma_D towards 1. The long run and the two transitions
  touching it are dropped and the alignment continues as two independent paths.
  Dropping the whole read instead would bias the other way, since reads crossing
  hard sequence carry most of the short gaps too.
- **I<->D adjacencies have no edge in this model.** They are counted and
  reported, never folded into a neighbour.

## 4. What was found

{c['reads']:,} alignments, {c['m_bases'] / 1e9:.1f} Gb of aligned columns.

```
{chr(10).join(ev)}
```

Two mismatch figures are given because they are measured differently. The NM
figure subtracts every indel base from the aligner's edit distance, over all
{c['reads']:,} reads. The base-by-base figure compares read to reference
directly over a {cal['reads']:,}-read subsample, and is the one to trust: it
cannot be confused by an edit distance that also counts something else. The two
differ by {100 * abs(c['mismatches'] / c['m_bases'] - mm_rate) / mm_rate:.0f}%.

Deletions open {p['delta_d'] / p['delta_i']:.1f}x more often than insertions,
and both are short: {100 * (1 - p['gamma_i']):.0f}% of insertions and
{100 * (1 - p['gamma_d']):.0f}% of deletions are a single base.

![Observed run lengths against the geometric law the model assumes.]({figs['runs']})

The bars leave the line at the tail. The model can only produce geometric runs,
so the excess of long gaps is a misspecification, not an estimation error — it
is what context-dependent gap parameters would absorb.

## 5. The estimates

With the path known, the maximum-likelihood estimate of a transition probability
is the ratio of counts — no optimiser, no EM:

```
delta_I = n(M->I) / [n(M->M) + n(M->I) + n(M->D)]      gamma_I = n(I->I) / [n(I->I) + n(I->M)]
delta_D = n(M->D) / [n(M->M) + n(M->I) + n(M->D)]      gamma_D = n(D->D) / [n(D->D) + n(D->M)]
```

```
{chr(10).join(est)}
```

**The jackknife is the honest error bar.** The binomial standard error treats
every aligned base as an independent draw, which they are not: gap rates vary
with context and reads in one tile share it. Deleting one tile at a time and
re-estimating gives a spread {jk['delta_i'] / d['se']['delta_i']:.1f}x wider,
over {jk['windows']:,} tiles.

![Each tile gives its own estimate; the histogram's width is what the jackknife measures.]({figs['tiles']})

Everything else in the model is fixed or derived rather than fitted. M->M is
1 - delta_I - delta_D - tau, B->I is delta_I, I->M is 1 - gamma_I - tau. The
termination probability tau is not estimated: read length is conditioned on, and
fitting tau to observed lengths would import the length model that conditioning
exists to remove. The per-base error probability eps comes from the quality
string, not from this fit — which section 6 shows to be a problem.

## 6. Emission calibration, measured but not applied

The same pass compares read to reference base by base and tallies by reported
quality. Reported qualities are pessimistic at every level:

```
{chr(10).join(table(qrows))}
```

![Every point sits below the diagonal: bases are better than they claim.]({figs['calib']})

`eps_from_phred` takes the reported value at face value, so a substitution
currently costs log(3(1-eps)/eps) = 10.3 nats where the data say 11.9 — a ~15%
under-charge wherever a locus is discriminated by substitutions rather than
indels, and a much larger one on low-quality bases, which the model presently
treats as nearly uninformative.

This is **not** applied, because the measurement is confounded in both
directions: the curated assembly was built from these same reads, so a read
error that reached the consensus reads as a match, deflating the observed rate;
residual assembly errors and haplotype mismapping read as read errors, inflating
it. A clean estimate wants a held-out assembly or variant-masked regions.

## 7. Context dependence

![Gap rate inside homopolymers against everywhere else.]({figs['context']})

Gaps are {hp / other:.1f}x more frequent inside homopolymers of 4 bp or more
({hp:.2e} per base) than elsewhere ({other:.2e} per base). One global gap rate
is therefore wrong in a measurable way, and wrong exactly where HiFi indel error
concentrates. This is the size of what context-dependent gap parameters would
buy, and the reason the seed locus's microsatellite verdict is the least
trustworthy magnitude the model produces.

## 8. What this cannot tell you

- **The path is minimap2's, not the model's posterior.** Its affine penalties
  (`map-hifi`: O 4,24 / E 2,1) decide which gaps exist at all, and its leftmost
  convention decides where a tied gap is placed. The estimator itself is checked
  by a round trip on sampled reads, where the path really is the model's, and
  recovers the generating parameters within 4 standard errors on all four.
- **The reference is not independent of the reads.** See section 6.
- **{f['max_gap']} bp as the gap ceiling is a judgement.** Sequencing indels and
  small variants overlap in size, and nothing in the data marks the boundary.
- **Fitting against the other assembly moves delta_I by 5%** — 60x the
  jackknife error. Which assembly is used as the reference, not the sample size,
  is what limits this estimate.
- **Contigs shorter than one tile are unrepresented.**
"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="docs/phmm_params_hg002.json")
    ap.add_argument("--bam", default="examples/alignments/hg002.curated.sorted.bam")
    ap.add_argument("--out", default="docs/phmm_params_report.pdf")
    ap.add_argument("--figdir", default="docs/figs")
    a = ap.parse_args(argv)

    with open(a.json) as fh:
        d = json.load(fh)
    os.makedirs(a.figdir, exist_ok=True)
    figs = {k: os.path.join(a.figdir, f"phmm_params_{k}.png")
            for k in ("runs", "calib", "context", "tiles")}

    per_tile = per_tile_counts(a.bam, os.path.join(a.figdir, "per_tile.json"))
    fig_run_lengths(d["counts"], d["params"], figs["runs"])
    fig_calibration(d["calibration"], figs["calib"])
    fig_context(d["calibration"], figs["context"])
    fig_between_tiles(per_tile, d["jackknife"], figs["tiles"])

    md_path = os.path.splitext(a.out)[0] + ".md"
    with open(md_path, "w") as fh:
        fh.write(build_markdown(d, per_tile, figs, a.bam))
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "md_to_pdf.py"),
                    md_path, a.out, "--title", "Pair-HMM parameters from HiFi alignments"],
                   check=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
