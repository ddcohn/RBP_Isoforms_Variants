import matplotlib.pyplot as plt

categories = [
    "COSMIC gene ID mapping",
    "Tissue expression",
    "Disease associations",
    "Core gene/transcript IDs",
    "RNA / CDS sequences",
    "GO terms",
    "STRING PPI",
]
values = [98.3, 97.7, 96.8, 95.8, 95.7, 86.1, 81.0]

# reverse so the highest value plots at the top of a horizontal bar chart
categories = categories[::-1]
values = values[::-1]

fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
bars = ax.barh(categories, values, color="#4472C4")  # standard Excel-blue

ax.set_xlim(0, 105)
ax.set_xlabel("% of rows populated")
ax.set_title("Data Coverage by Category")

for bar, val in zip(bars, values):
    ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2,
             f"{val:.1f}%", va="center", fontsize=10)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("/u/home/d/ddcohn/_claude_plain_chart.png")
print("saved")
