# Auditing the lab's existing ClinVar variant-stats table

Not this project's own pipeline — a diagnostic pass over a labmate's
(`tchhabri`) already-built ClinVar variant-classification table
(`Isoform_Post_Merge_PSLab_OpenTargets_withVariantStats.csv`), prompted by
a labmate suggesting COSMIC as an additional data source and a general
"does this table's data hold up" concern. See `../FINDINGS.md` for the
full narrative; this directory has just the two scripts that produced the
numbers.

## What we found

1. **`variant_stats_summary_v1.py`** — naive per-row sum of
   `Total_Pathogenic` / `Total_VUS` / `Total_Benign` across the whole
   65,744-row table. Result: **596.7M pathogenic + 497.2M VUS + 596.5M
   benign — 1.69 billion classified variant instances**, absurd for a
   table covering ~17,582 distinct proteins (ClinVar's *entire* database
   is only ~3-4M variants total). Traced the cause: the table has the same
   duplicate-row problem as the original broken isoform table (`BRCA1`
   appears 369 times), and this script summed across every duplicate row.

2. **`variant_stats_summary_v2.py`** — dedupes to one row per protein
   first (17,582 distinct proteins), and confirms the duplicate rows for
   any given protein carry byte-identical stats (0/17,582 disagree) — so
   the row duplication itself isn't corrupting the per-protein numbers,
   it's a separate problem stacked on top. Even after deduping, the real
   total is **32.7M classified variant instances across 17,582 proteins**
   — still roughly 10x too large. Root cause (documented but not yet
   fixed): the pipeline's own documentation states "the same variant
   appears multiple times when annotated on multiple isoforms — dedupe by
   `VariationID` before counting," and that dedup step appears to be
   missing from whatever script rolled per-variant ClinVar records up into
   these per-protein totals.

## Where the row duplication actually comes from

Traced via `merge_stats_to_isoform.py` (in `tchhabri`'s own directory, not
copied here) to a **third, separate isoform-table lineage**
(`/u/project/kappel/RBP/Isoform_Table/Old/Isoform_Post_Merge_PSLab_OpenTargets.csv`,
owned by `fraza`, dated April) — distinct from both this project's clean
rebuild and the originally-diagnosed broken table. `Dominant_Isoform` is
inherited from that file unchanged; nothing in `tchhabri`'s pipeline
computes it. See `../FINDINGS.md`.
