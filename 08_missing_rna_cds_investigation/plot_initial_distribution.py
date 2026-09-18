import matplotlib.pyplot as plt

categories = [
    "Antisense/lncRNA/\nuncharacterized locus",
    "Other/unclassified\n(no exact match found)",
    "Pseudogene",
    "No ensembl_gene\nat all",
    "Immunoglobulin\ngene segment",
    "Mitochondrial-encoded\nmicropeptide",
    "TCR gene segment",
]
counts = [318, 242, 142, 140, 16, 12, 6]

categories = categories[::-1]
counts = counts[::-1]

fig, ax = plt.subplots(figsize=(8.5, 5), dpi=150)
bars = ax.barh(categories, counts, color="#4472C4")

ax.set_xlabel("Number of proteins (of 876 missing RNA/CDS)")
ax.set_title("Why 876 Proteins Are Missing RNA/CDS Sequence")

for bar, val in zip(bars, counts):
    pct = 100 * val / 876
    ax.text(bar.get_width() + 4, bar.get_y() + bar.get_height() / 2,
             f"{val} ({pct:.1f}%)", va="center", fontsize=9.5)

ax.set_xlim(0, 360)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("/u/home/d/ddcohn/_claude_missing_cds_chart.png")
print("saved")
