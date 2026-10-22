import json
import statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500"]
INK = "#023047"
MUTED = "#6b7280"

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.edgecolor": "#d1d5db",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

with open("/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/satmut_vs_clinvar_cosmic.json") as f:
    data = json.load(f)

buckets = data["buckets"]
n_total = data["n_total"]
counts = data["counts_matched"]

BUCKET_ORDER = ["background", "cosmic", "clinvar_benign", "clinvar_pathogenic"]
BUCKET_LABELS = {
    "background": "All mutations\n(background)",
    "cosmic": "COSMIC\n(no sig. label)",
    "clinvar_benign": "ClinVar:\nBenign",
    "clinvar_pathogenic": "ClinVar:\nPathogenic",
}
BUCKET_COLORS = {
    "background": MUTED,
    "cosmic": PALETTE[3],
    "clinvar_benign": PALETTE[1],
    "clinvar_pathogenic": PALETTE[4],
}

METRICS = [
    ("tango", "TANGO |Δ AGG| vs WT"),
    ("protgps", "protGPS total |Δ score| vs WT (12 compartments)"),
    ("deeploc", "DeepLoc total |Δ score| vs WT (14 compartments)"),
]

fig, axes = plt.subplots(1, 3, figsize=(17, 6.8))
fig.subplots_adjust(wspace=0.32, top=0.76, bottom=0.16, left=0.06, right=0.98)

means_by_metric = {}
for ax, (metric, mlabel) in zip(axes, METRICS):
    box_data = [buckets[b][metric] for b in BUCKET_ORDER]
    nonzero = [v for arr in box_data for v in arr if v > 0]
    linthresh = max(statistics.median(nonzero) * 0.1, 1e-4) if nonzero else 1e-4

    bp = ax.boxplot(
        box_data, tick_labels=[BUCKET_LABELS[b] for b in BUCKET_ORDER],
        showfliers=False, patch_artist=True, widths=0.55, medianprops={"color": INK, "linewidth": 1.5},
    )
    for patch, b in zip(bp["boxes"], BUCKET_ORDER):
        patch.set_facecolor(BUCKET_COLORS[b])
        patch.set_alpha(0.85)
        patch.set_edgecolor(INK)
    for element in ["whiskers", "caps"]:
        for line in bp[element]:
            line.set_color(MUTED)

    ax.set_yscale("symlog", linthresh=linthresh)
    ax.set_title(mlabel, fontsize=10.5, pad=10)
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=9)

    means_by_metric[metric] = {b: (statistics.mean(v) if v else 0) for b, v in zip(BUCKET_ORDER, box_data)}

SHORT_BUCKET_NAME = {"background": "All", "cosmic": "COSMIC", "clinvar_benign": "ClinVar-Benign", "clinvar_pathogenic": "ClinVar-Pathogenic"}
SHORT_TOOL_NAME = {"tango": "TANGO", "protgps": "protGPS", "deeploc": "DeepLoc"}
caption_parts = []
for metric, _ in METRICS:
    vals = ", ".join(f"{SHORT_BUCKET_NAME[b]} {means_by_metric[metric][b]:.3g}" for b in BUCKET_ORDER)
    caption_parts.append(f"{SHORT_TOOL_NAME[metric]} mean: {vals}")
caption = "   |   ".join(caption_parts)

fig.suptitle(
    f"Saturation-scan prediction change vs WT: ClinVar + COSMIC real variants vs background\n"
    f"ClinVar Pathogenic n={counts['clinvar_pathogenic']:,}, Benign n={counts['clinvar_benign']:,} (curated labels)  |  "
    f"COSMIC n={counts['cosmic']:,} (all observed cancer missense matches -- no pathogenicity tier available in this data, NOT curated)\n"
    f"(symmetric-log y-axis: distributions are heavily right-skewed, most mutations change very little)",
    fontsize=11, x=0.01, ha="left",
)
fig.text(0.06, 0.03, caption, fontsize=7.5, color=MUTED)
fig.savefig("/Users/danielcohn/Desktop/satmut_vs_clinvar_cosmic_comparison.png", dpi=150)
plt.close(fig)

print("Wrote /Users/danielcohn/Desktop/satmut_vs_clinvar_cosmic_comparison.png")
for metric, mlabel in METRICS:
    print(f"\n{mlabel}:")
    for b in BUCKET_ORDER:
        vals = buckets[b][metric]
        if vals:
            print(f"  {b}: n={len(vals)}, mean={statistics.mean(vals):.4f}, median={statistics.median(vals):.4f}")
