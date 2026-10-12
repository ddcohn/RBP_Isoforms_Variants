"""Step 10: example figure of per-base phyloP across one RBP's coding sequence, coloured by codon
position (reads results/per_base/rbp_cds_per_base_scores.parquet from step 9).

Usage: python src/10_phylop_example_figure.py [GENE]   (default HNRNPA1)
Output: results/figures/<GENE>_phylop_per_base.png
Run from the 13_evolutionary_conservation directory.
"""
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PAL = dict(light="#8ecae6", mid="#219ebc", dark="#023047", gold="#ffb703", orange="#fb8500")
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                     "axes.edgecolor": PAL["dark"], "axes.labelcolor": PAL["dark"], "text.color": PAL["dark"],
                     "xtick.color": PAL["dark"], "ytick.color": PAL["dark"], "savefig.dpi": 300})


def main(gene="HNRNPA1"):
    b = pd.read_parquet("results/per_base/rbp_cds_per_base_scores.parquet",
                        columns=["GeneSymbol", "cds_pos", "codon_position", "exon_rank", "phyloP100way"])
    g = b[b["GeneSymbol"] == gene].sort_values("cds_pos")
    assert len(g), f"{gene} not in the per-base table"
    pos, cp, v = g["cds_pos"].to_numpy() - 1, g["codon_position"].to_numpy(), g["phyloP100way"].to_numpy()
    edges = pos[np.r_[False, np.diff(g["exon_rank"].to_numpy()) != 0]]
    n = len(pos)
    m12, m3, mg = np.nanmean(v[cp < 3]), np.nanmean(v[cp == 3]), np.nanmean(v)

    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    fig.subplots_adjust(left=0.10, right=0.80, top=0.86, bottom=0.30)
    for e in edges:
        ax.axvline(e, color=PAL["light"], lw=0.6, zorder=0)
    ax.axhline(0, color=PAL["dark"], lw=0.5, zorder=1)
    ax.scatter(pos[cp < 3], v[cp < 3], s=3.5, color=PAL["mid"], alpha=0.75, linewidths=0, zorder=2)
    ax.scatter(pos[cp == 3], v[cp == 3], s=3.5, color=PAL["orange"], alpha=0.85, linewidths=0, zorder=3)
    for y, c, t, ty, ls in [(m12, PAL["mid"], f"positions 1 and 2: {m12:.1f}", max(m12 + 1.1, 7.0) if m12 > 5 else m12, "-"),
                            (mg, PAL["dark"], f"whole gene: {mg:.1f}", mg, "--"),
                            (m3, PAL["orange"], f"position 3: {m3:.1f}", m3, "-")]:
        ax.hlines(y, 0, n, color=c, lw=1.1, ls=ls, zorder=4)
        ax.text(n * 1.015, ty, t, color=c, va="center", fontsize=7)
    ax.set_xlim(-10, n + 10)
    lo = np.floor(np.nanmin(v)) - 0.5  # show every base: no points are clipped by the axis
    ax.set_ylim(lo, 11)
    ax.set_xlabel(f"Position along the {gene} coding sequence (nt); light lines mark exon boundaries")
    ax.set_ylabel("phyloP per base\n(100 vertebrates)")
    ax.set_title("Codon positions 1 and 2 are strongly conserved; position 3 is near neutral"
                 if m12 - m3 > 2 else f"Per-base phyloP across {gene}", loc="left", pad=6)
    ax.set_yticks([t for t in (-10, -5, 0, 2, 5, 10) if t >= lo])
    fig.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", ms=4, color=PAL["mid"], label="Codon positions 1 and 2"),
                        mpl.lines.Line2D([], [], marker="o", ls="", ms=4, color=PAL["orange"], label="Codon position 3")],
               loc="lower left", bbox_to_anchor=(0.10, 0.0), ncol=2, frameon=False, handletextpad=0.2)
    out = Path("results/figures") / f"{gene}_phylop_per_base.png"
    fig.savefig(out)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    # tick labels for ticks outside the view are never drawn; leave them out of the check
    hidden = {id(tk.label1) for axis in (ax.xaxis, ax.yaxis)
              for tk in axis.get_major_ticks()
              if not min(axis.get_view_interval()) <= tk.get_loc() <= max(axis.get_view_interval())}
    boxes = [(t.get_text(), t.get_window_extent(rend)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and id(t) not in hidden]
    off = [t for t, b in boxes if not (b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.x1 and b.y1 <= fig.bbox.y1)]
    assert not off, f"text outside the canvas: {off}"
    print(out, f"n={n} pos1+2={m12:.2f} pos3={m3:.2f} gene={mg:.2f}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
