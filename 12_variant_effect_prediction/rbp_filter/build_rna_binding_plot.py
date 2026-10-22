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

N_TOTAL = 2335
N_GENERAL = 1286

MF_ONLY = {
    "mRNA": 253, "tRNA": 63, "rRNA": 60, "snRNA": 41, "miRNA": 39,
    "snoRNA": 23, "lncRNA": 17, "siRNA": 7, "piRNA": 6,
}
RECOVERED = {
    "mRNA": 228, "rRNA": 121, "tRNA": 32, "miRNA": 15, "snRNA": 13,
    "snoRNA": 11, "lncRNA": 3, "piRNA": 5, "siRNA": 4,
}
N_STILL_GAP = 603

order = sorted(MF_ONLY, key=lambda k: -(MF_ONLY[k] + RECOVERED[k]))
labels = order + ["RNA binding\n(still no class)"]
mf_vals = [MF_ONLY[k] for k in order] + [0]
rec_vals = [RECOVERED[k] for k in order] + [N_STILL_GAP]

fig, ax = plt.subplots(figsize=(11, 6.5))
x = range(len(labels))
bars1 = ax.bar(x, mf_vals, color=PALETTE[1], width=0.6, label="From GO Molecular Function (original)")
bars2 = ax.bar(x, rec_vals, bottom=mf_vals, color=PALETTE[3], width=0.6, label="Recovered via GO Biological Process + HGNC family")
colors_gap = [PALETTE[1]] * len(order) + [MUTED]
bars2[-1].set_color(MUTED)

ax.set_xticks(list(x))
ax.set_xticklabels(labels)
ax.set_title(
    f"RNA-binding GO annotation coverage across the WT RBP list\n"
    f"{N_GENERAL:,} of {N_TOTAL:,} (55.1%) have general 'RNA binding'; specific-class coverage now {816:,} ({816/N_TOTAL*100:.1f}%), up from {450:,} ({450/N_TOTAL*100:.1f}%)",
    fontsize=10.5, loc="left",
)
ax.set_ylabel("Protein count")
ax.grid(axis="y", color="#e5e7eb", linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=9, loc="upper left")

for i, (mf, rec) in enumerate(zip(mf_vals, rec_vals)):
    total = mf + rec
    if total > 0:
        ax.text(i, total + 10, f"{total:,}", ha="center", fontsize=9, color=INK)

ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
fig.text(0.01, 0.01, "Specific-class bars are not mutually exclusive. Excluded one real false positive during curation: HGNC family\n'MicroRNA protein coding host genes' (159 proteins) means the gene's own intron hosts a miRNA gene -- not that the protein binds miRNA.", fontsize=7.5, color=MUTED)
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig("/Users/danielcohn/Desktop/wt_rbp_rna_binding_go_breakdown.png", dpi=150)
print("Wrote /Users/danielcohn/Desktop/wt_rbp_rna_binding_go_breakdown.png")
