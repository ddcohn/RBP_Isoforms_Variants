import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

PALETTE = ["#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500"]
INK = "#023047"
MUTED = "#6b7280"

with open("/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/satmut_plot_data.json") as f:
    data = json.load(f)

n_total = data["n_total"]

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

# --- Plot 1: TANGO AGG score distribution ---
fig, ax = plt.subplots(figsize=(8, 5))
agg = data["tango_agg_sample"]
ax.hist(agg, bins=60, color=PALETTE[1], edgecolor="white", linewidth=0.3)
ax.set_title(f"TANGO aggregation propensity (AGG) across all point mutants\nn = {len(agg):,} sampled of {n_total:,} total", fontsize=10.5, loc="left")
ax.set_xlabel("AGG score (whole-sequence sum)")
ax.set_ylabel("Mutant count")
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.grid(axis="y", color="#e5e7eb", linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig("/Users/danielcohn/Desktop/satmut_tango_agg_distribution.png", dpi=150)
plt.close(fig)

# --- Plot 2: protGPS mean score per compartment ---
means = data["protgps_means"]
na_counts = data["protgps_na_counts"]
na_n = list(na_counts.values())[0]
items = sorted(means.items(), key=lambda kv: kv[1], reverse=True)
labels = [k.replace("_Score", "").replace("_", " ") for k, v in items]
values = [v for k, v in items]

fig, ax = plt.subplots(figsize=(9.5, 5.5))
bars = ax.barh(labels, values, color=PALETTE[1], height=0.65)
ax.invert_yaxis()
ax.set_title(
    f"protGPS: mean predicted condensate-partitioning probability per mutant\n"
    f"n = {n_total - na_n:,} scored ({na_n:,} excluded: sequence too long for the model)",
    fontsize=10.5, loc="left",
)
ax.set_xlabel("Mean predicted probability")
ax.grid(axis="x", color="#e5e7eb", linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
for bar, v in zip(bars, values):
    ax.text(v + 0.004, bar.get_y() + bar.get_height() / 2, f"{v:.3f}", va="center", fontsize=9, color=INK)
fig.tight_layout()
fig.savefig("/Users/danielcohn/Desktop/satmut_protgps_mean_scores.png", dpi=150)
plt.close(fig)

# --- Plot 3: DeepLoc localization frequency ---
loc_counts = data["deeploc_loc_counts"]
items = sorted(loc_counts.items(), key=lambda kv: kv[1], reverse=True)
labels = [k for k, v in items]
values = [v for k, v in items]
pcts = [100.0 * v / n_total for v in values]

fig, ax = plt.subplots(figsize=(9.5, 5.5))
bars = ax.barh(labels, values, color=PALETTE[1], height=0.65)
ax.invert_yaxis()
ax.set_title(
    f"DeepLoc2.1: predicted subcellular localization frequency\n"
    f"n = {n_total:,} mutants (can have >1 predicted localization each)",
    fontsize=10.5, loc="left",
)
ax.set_xlabel("Mutant count")
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.grid(axis="x", color="#e5e7eb", linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
for bar, v, p in zip(bars, values, pcts):
    ax.text(v + n_total * 0.005, bar.get_y() + bar.get_height() / 2, f"{p:.1f}%", va="center", fontsize=9, color=INK)
fig.tight_layout()
fig.savefig("/Users/danielcohn/Desktop/satmut_deeploc_localization_frequency.png", dpi=150)
plt.close(fig)

print("Wrote 3 PNGs to /Users/danielcohn/Desktop/")
print("protGPS means (sorted):", items if False else dict(sorted(means.items(), key=lambda kv: kv[1], reverse=True)))
print("DeepLoc counts (sorted):", dict(sorted(loc_counts.items(), key=lambda kv: kv[1], reverse=True)))
