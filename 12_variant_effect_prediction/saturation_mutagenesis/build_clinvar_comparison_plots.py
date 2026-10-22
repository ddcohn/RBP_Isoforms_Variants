import json
import statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

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

with open("/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/satmut_vs_clinvar.json") as f:
    data = json.load(f)

buckets = data["buckets"]
n_total = data["n_total"]
n_matched = data["n_clinvar_matched"]

BUCKET_ORDER = ["background", "uncertain", "benign", "pathogenic"]
BUCKET_LABELS = {
    "background": "All mutations\n(background)",
    "uncertain": "ClinVar:\nUncertain sig.",
    "benign": "ClinVar:\nBenign",
    "pathogenic": "ClinVar:\nPathogenic",
}
BUCKET_COLORS = {
    "background": MUTED,
    "uncertain": PALETTE[0],
    "benign": PALETTE[1],
    "pathogenic": PALETTE[4],
}

METRICS = [
    ("tango", "TANGO |Δ AGG| vs WT"),
    ("protgps", "protGPS total |Δ score| vs WT (12 compartments)"),
    ("deeploc", "DeepLoc total |Δ score| vs WT (14 compartments)"),
]

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
fig.subplots_adjust(wspace=0.3, top=0.82, bottom=0.12, left=0.06, right=0.98)

for ax, (metric, mlabel) in zip(axes, METRICS):
    box_data = [buckets[b][metric] for b in BUCKET_ORDER]
    bp = ax.boxplot(
        box_data, labels=[BUCKET_LABELS[b] for b in BUCKET_ORDER],
        showfliers=False, patch_artist=True, widths=0.55, medianprops={"color": INK, "linewidth": 1.5},
    )
    for patch, b in zip(bp["boxes"], BUCKET_ORDER):
        patch.set_facecolor(BUCKET_COLORS[b])
        patch.set_alpha(0.85)
        patch.set_edgecolor(INK)
    for element in ["whiskers", "caps"]:
        for line in bp[element]:
            line.set_color(MUTED)
    ax.set_title(mlabel, fontsize=10.5)
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=9)

    means = [statistics.mean(v) if v else 0 for v in box_data]
    for i, m in enumerate(means, start=1):
        ax.text(i, ax.get_ylim()[1] * 0.97, f"mean\n{m:.3g}", ha="center", va="top", fontsize=7.5, color=INK)

fig.suptitle(
    f"Saturation-scan prediction change vs WT: real ClinVar missense variants vs background\n"
    f"{n_matched:,} of {n_total:,} satmut mutations ({100*n_matched/n_total:.1f}%) correspond to a real observed ClinVar missense variant",
    fontsize=12, x=0.01, ha="left",
)
fig.savefig("/Users/danielcohn/Desktop/satmut_vs_clinvar_comparison.png", dpi=150)
plt.close(fig)

print("Wrote /Users/danielcohn/Desktop/satmut_vs_clinvar_comparison.png")
print()
for metric, mlabel in METRICS:
    print(f"\n{mlabel}:")
    for b in BUCKET_ORDER:
        vals = buckets[b][metric]
        if vals:
            print(f"  {b}: n={len(vals)}, mean={statistics.mean(vals):.4f}, median={statistics.median(vals):.4f}")
