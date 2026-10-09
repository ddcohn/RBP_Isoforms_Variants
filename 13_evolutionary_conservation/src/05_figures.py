"""Step 5: summary figure of RBP CDS conservation.  Run from the
13_evolutionary_conservation directory.  Output: results/figures/rbp_conservation_summary.png/.pdf
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("results")
PAL = {"light": "#8ecae6", "mid": "#219ebc", "dark": "#023047", "gold": "#ffb703", "orange": "#fb8500"}
LADDER = ["ptroglodytes", "mmulatta", "mmusculus", "rnorvegicus", "ggallus",
          "xtropicalis", "drerio", "dmelanogaster", "celegans", "scerevisiae"]
LABEL = {"ptroglodytes": "Chimpanzee", "mmulatta": "Rhesus macaque", "mmusculus": "Mouse",
         "rnorvegicus": "Rat", "ggallus": "Chicken", "xtropicalis": "Xenopus",
         "drerio": "Zebrafish", "dmelanogaster": "Fly", "celegans": "Worm",
         "scerevisiae": "Yeast"}
CLASSES = ["Spliceosomal", "Ribosomal", "Other RNP machine", "Other RBPs"]
TICK = {"Spliceosomal": "Spliceosomal", "Ribosomal": "Ribosomal",
        "Other RNP machine": "Other RNP\nmachine", "Other RBPs": "Other RBPs"}
CCOL = {"Spliceosomal": PAL["orange"], "Ribosomal": PAL["mid"],
        "Other RNP machine": PAL["gold"], "Other RBPs": PAL["light"]}

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.edgecolor": PAL["dark"],
    "axes.labelcolor": PAL["dark"], "text.color": PAL["dark"], "xtick.color": PAL["dark"],
    "ytick.color": PAL["dark"], "figure.dpi": 150, "savefig.dpi": 300, "pdf.fonttype": 42,
})


def classify(r):
    # precedence for RBPs flagged in more than one RNP-machine column
    if r["Spliceosomal_protein"] == "Y":
        return "Spliceosomal"
    if r["Ribosomal_protein"] == "Y":
        return "Ribosomal"
    if r["Other_RNP_machine"] == "Y":
        return "Other RNP machine"
    return "Other RBPs"


def main():
    t = pd.read_csv(OUT / "tables/rbp_conservation_scores.tsv", sep="\t")
    t = t[t["match_type"] != "none"].copy()
    t["cls"] = t.apply(classify, axis=1)
    n = len(t)
    counts = t["cls"].value_counts().to_dict()

    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 1.1, 1.0], wspace=0.5,
                          left=0.065, right=0.985, top=0.78, bottom=0.30)
    axa, axb, axc = (fig.add_subplot(g) for g in gs)

    # a: distribution of CDS phyloP by class (jitter + median tick)
    rng = np.random.default_rng(0)
    for i, c in enumerate(CLASSES):
        v = t.loc[t["cls"] == c, "phyloP100way_mean"].to_numpy()
        axa.scatter(i + rng.uniform(-0.28, 0.28, v.size), v, s=2.5, color=CCOL[c],
                    alpha=0.55, linewidths=0, rasterized=True)
        axa.hlines(np.median(v), i - 0.36, i + 0.36, color=PAL["dark"], lw=1.4, zorder=3)
    axa.set_xticks(range(len(CLASSES)))
    axa.set_xticklabels([f"{TICK[c]}\nn={counts.get(c, 0)}" for c in CLASSES])
    axa.tick_params(axis="x", labelsize=6)
    axa.set_ylabel("Mean CDS phyloP, 100 vertebrates")
    axa.set_title("Spliceosomal RBPs are most conserved", loc="left", pad=6)
    axa.margins(x=0.04, y=0.04)

    # b: vertebrate vs primate conservation, label extremes
    for c in CLASSES[::-1]:
        s = t[t["cls"] == c]
        axb.scatter(s["phyloP100way_mean"], s["phyloP30way_mean"], s=3,
                    color=CCOL[c], alpha=0.6, linewidths=0, label=c, rasterized=True)
    r = np.corrcoef(t["phyloP100way_mean"], t["phyloP30way_mean"])[0, 1]
    axb.set_xlabel("Mean CDS phyloP, 100 vertebrates")
    axb.set_ylabel("Mean CDS phyloP, 30-way (primates)")
    axb.set_title("Vertebrate and primate scores agree", loc="left", pad=6)
    axb.text(0.04, 0.96, f"Pearson r = {r:.2f}", transform=axb.transAxes, va="top", fontsize=7)
    lo = t.nsmallest(2, "phyloP100way_mean").sort_values("phyloP30way_mean", ascending=False)
    for (_, row), dy in zip(lo.iterrows(), (9, -11)):  # stagger the two neighbours
        axb.annotate(row["GeneSymbol"], (row["phyloP100way_mean"], row["phyloP30way_mean"]),
                     xytext=(14, dy), textcoords="offset points", fontsize=6, style="italic",
                     va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=PAL["dark"]))
    hi = t.nlargest(1, "phyloP100way_mean").iloc[0]
    axb.annotate(hi["GeneSymbol"], (hi["phyloP100way_mean"], hi["phyloP30way_mean"]),
                 xytext=(-14, -16), textcoords="offset points", fontsize=6, style="italic",
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=PAL["dark"]))
    axb.margins(0.06)

    # c: share of RBPs with a Compara ortholog, ordered by divergence from human
    frac = [t[f"{sp}_orthology_type"].notna().mean() * 100 for sp in LADDER]
    y = np.arange(len(LADDER))[::-1]
    axc.hlines(y, 0, frac, color=PAL["light"], lw=1.2)
    axc.scatter(frac, y, s=18, color=PAL["mid"], zorder=3)
    for yi, f in zip(y, frac):
        axc.text(f + 2.5, yi, f"{f:.0f}", ha="left", va="center", fontsize=6, color=PAL["dark"])
    axc.set_yticks(y)
    axc.set_yticklabels([LABEL[s] for s in LADDER])
    axc.set_xlim(0, 112)
    axc.set_xticks([0, 25, 50, 75, 100])
    axc.set_xlabel("RBPs with an ortholog (%)")
    axc.set_title("Orthology thins out\nbeyond vertebrates", loc="left", pad=6)
    axc.margins(y=0.05)

    # one shared class legend, below panels a-b
    handles = [mpl.lines.Line2D([], [], marker="o", ls="", ms=4, color=CCOL[c], label=c)
               for c in CLASSES]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.07, 0.0), ncol=4,
               frameon=False, handletextpad=0.2, columnspacing=1.2)
    for ax, letter in zip((axa, axb, axc), "abc"):
        ax.text(-0.02, 1.24, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
                va="bottom", ha="right")
    (OUT / "figures").mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"figures/rbp_conservation_summary.{ext}")
    print("n RBPs:", n, counts, "r =", round(r, 3))
    print("ortholog %:", dict(zip(LADDER, np.round(frac, 1))))
    print("median phyloP100 by class:",
          t.groupby("cls")["phyloP100way_mean"].median().round(2).to_dict())


if __name__ == "__main__":
    main()
