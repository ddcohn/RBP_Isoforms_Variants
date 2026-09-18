import matplotlib.pyplot as plt

categories = [
    "Pseudogene",
    "Antisense/lncRNA/\nuncharacterized locus",
    "No ensembl_gene\nat all",
    "Ensembl ID not found\n/ retired",
    "Mitochondrial-encoded\nmicropeptide",
    "Immunoglobulin\ngene segment",
    "TCR gene segment",
]
counts = [311, 267, 278, 9, 8, 2, 1]

categories = categories[::-1]
counts = counts[::-1]

fig, ax = plt.subplots(figsize=(8.5, 5), dpi=150)
bars = ax.barh(categories, counts, color="#4472C4")

ax.set_xlabel("Number of proteins (of 876 missing RNA/CDS)")
ax.set_title("Why 876 Proteins Are Missing RNA/CDS Sequence\n(verified against Ensembl biotype)")

for bar, val in zip(bars, counts):
    pct = 100 * val / 876
    ax.text(bar.get_width() + 4, bar.get_y() + bar.get_height() / 2,
             f"{val} ({pct:.1f}%)", va="center", fontsize=9.5)

ax.set_xlim(0, 340)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("/u/home/d/ddcohn/_claude_corrected_chart.png")
print("saved")
