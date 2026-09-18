# Column-by-column stats and inventory

Generic tooling used repeatedly throughout this project whenever "what
does the current table actually look like" needed a real answer instead
of a guess.

- **`table_stats_v1.py` → `v2.py` → `v3.py`** — iterative versions of a
  per-column summary script (% populated, distinct value counts, min/max/
  mean/median for numeric or length-like columns, categorical breakdowns
  for enum-like columns). Each version added a data type the previous one
  handled naively — v1 was a first pass on core columns, v2 added
  proper numeric-vs-zero handling (a column full of `"0.0"` strings was
  initially miscounted as "populated" when it should count as "no IDR"),
  v3 added semicolon-separated list stats for the PPI/GO/OpenTargets-style
  columns.

- **`full_inventory.py`** — runs the v3-style stats across *two* tables at
  once (the isoform table and the COSMIC gene summary) and writes one
  combined TSV, for a single Excel-pasteable overview instead of two
  separate ones.

- **`plot_coverage_summary.py`** — takes the "interesting" (non-100%)
  columns from the inventory above, groups them into logical categories
  (Core IDs, RNA/CDS, PPI, GO, tissue, disease associations, COSMIC gene
  mapping), and renders a plain matplotlib bar chart of coverage by
  category. Written deliberately in plain matplotlib defaults (no custom
  styling/theming) per explicit request — a custom-styled HTML/CSS chart
  was built first and rejected as looking too obviously AI-generated.

`submit_inventory.sh` is the UGE job script — the isoform table is large
enough (500MB+) that computing per-column stats needs to run off the
login node.
