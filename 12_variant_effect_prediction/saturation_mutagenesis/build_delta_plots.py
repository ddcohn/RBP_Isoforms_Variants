import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap

SCRATCH = "/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad"
DESKTOP = "/Users/danielcohn/Desktop"

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

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

METRICS = [
    ("mean_tango_abs_delta_agg", "TANGO |Δ AGG|", PALETTE[1]),
    ("mean_protgps_total_abs_delta", "protGPS total |Δ score| (12 compartments)", PALETTE[3]),
    ("mean_deeploc_total_abs_delta", "DeepLoc total |Δ score| (14 compartments)", PALETTE[4]),
]


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# --- Plot 1 & 2: change by wt_aa / mut_aa ---
for fname, label_col, title_prefix, outname in [
    ("satmut_change_by_wt_aa.tsv", "wt_aa", "Original residue (WT)", "satmut_change_by_wt_aa.png"),
    ("satmut_change_by_mut_aa.tsv", "mut_aa", "Substituted residue (mutant)", "satmut_change_by_mut_aa.png"),
]:
    rows = {r[label_col]: r for r in read_tsv(f"{SCRATCH}/{fname}")}
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=False)
    for ax, (col, mlabel, color) in zip(axes, METRICS):
        vals = [float(rows[aa][col]) for aa in AA_ORDER]
        order = sorted(range(20), key=lambda i: -vals[i])
        labels_sorted = [AA_ORDER[i] for i in order]
        vals_sorted = [vals[i] for i in order]
        ax.bar(labels_sorted, vals_sorted, color=color, width=0.7)
        ax.set_title(mlabel, fontsize=10.5)
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelsize=9)
    fig.suptitle(f"{title_prefix}: mean prediction change vs WT, by amino acid\nn = 31,525,484 mutations total", fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{DESKTOP}/{outname}", dpi=150)
    plt.close(fig)

# --- Plot 3: substitution pair heatmap (protGPS) ---
pair_rows = read_tsv(f"{SCRATCH}/satmut_change_by_substitution_pair.tsv")
matrix = {(r["wt_aa"], r["mut_aa"]): float(r["mean_protgps_total_abs_delta"]) for r in pair_rows}
grid = [[matrix.get((wt, mut), float("nan")) for mut in AA_ORDER] for wt in AA_ORDER]

cmap = LinearSegmentedColormap.from_list("seq", ["#ffffff", PALETTE[0], PALETTE[1], PALETTE[2]])
fig, ax = plt.subplots(figsize=(9.5, 8.5))
im = ax.imshow(grid, cmap=cmap, aspect="auto")
ax.set_xticks(range(20))
ax.set_xticklabels(AA_ORDER, fontsize=9)
ax.set_yticks(range(20))
ax.set_yticklabels(AA_ORDER, fontsize=9)
ax.set_xlabel("Substituted residue (mutant)")
ax.set_ylabel("Original residue (WT)")
ax.set_title("protGPS total |Δ score| by substitution pair\n(white/light = little change, dark = large change)", fontsize=11, loc="left")
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xticks([x - 0.5 for x in range(21)], minor=True)
ax.set_yticks([y - 0.5 for y in range(21)], minor=True)
ax.grid(which="minor", color="white", linewidth=1.5)
ax.tick_params(which="minor", length=0)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Mean total |Δ score|", fontsize=9)
fig.tight_layout()
fig.savefig(f"{DESKTOP}/satmut_substitution_heatmap_protgps.png", dpi=150)
plt.close(fig)

# --- Plot 4: top susceptible proteins, one panel per metric ---
protein_rows = read_tsv(f"{SCRATCH}/satmut_change_by_protein.tsv")
protein_rows = [r for r in protein_rows if int(r["n_mutations"]) >= 1000]  # drop tiny proteins (few mutants -> noisy mean)

def short_label(raw):
    parts = raw.split(";")
    if len(parts) == 1:
        return parts[0]
    return f"{parts[0]} (+{len(parts) - 1} alias)"


fig, axes = plt.subplots(1, 3, figsize=(19, 7.5))
fig.subplots_adjust(wspace=0.55, top=0.84, bottom=0.08, left=0.08, right=0.97)
for ax, (col, mlabel, color) in zip(axes, METRICS):
    valid_rows = [r for r in protein_rows if r[col].strip()]
    ranked = sorted(valid_rows, key=lambda r: -float(r[col]))[:20]
    labels = [short_label(r["gene_symbol"] or r["uniprot_accession"]) for r in ranked][::-1]
    vals = [float(r[col]) for r in ranked][::-1]
    ax.barh(labels, vals, color=color, height=0.65)
    ax.set_title(mlabel, fontsize=10, wrap=True)
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4))

fig.suptitle("Top 20 proteins most susceptible to mutation-driven prediction change\n(proteins with >=1,000 scored mutants only)", fontsize=12, x=0.01, ha="left")
fig.savefig(f"{DESKTOP}/satmut_top_susceptible_proteins.png", dpi=150)
plt.close(fig)

print("Wrote 4 PNGs to Desktop")

# print some real numbers for the report
for col, mlabel, _ in METRICS:
    valid_rows = [r for r in protein_rows if r[col].strip()]
    ranked = sorted(valid_rows, key=lambda r: -float(r[col]))[:5]
    print(f"\nTop 5 by {mlabel}:")
    for r in ranked:
        print(f"  {r['gene_symbol']}: {float(r[col]):.4f} (n={r['n_mutations']})")

wt_rows = {r["wt_aa"]: r for r in read_tsv(f"{SCRATCH}/satmut_change_by_wt_aa.tsv")}
mut_rows = {r["mut_aa"]: r for r in read_tsv(f"{SCRATCH}/satmut_change_by_mut_aa.tsv")}
for col, mlabel, _ in METRICS:
    top_wt = max(wt_rows.items(), key=lambda kv: float(kv[1][col]))
    top_mut = max(mut_rows.items(), key=lambda kv: float(kv[1][col]))
    print(f"\n{mlabel}: most volatile WT residue = {top_wt[0]} ({float(top_wt[1][col]):.4f}), most disruptive substituted residue = {top_mut[0]} ({float(top_mut[1][col]):.4f})")
