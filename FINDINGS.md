# Findings that cut across multiple directories

Things discovered while doing this work that aren't specific to any one
pipeline stage, but matter for how the whole project (and the wider lab's
related files) should be understood going forward.

## The isoform table is not RBP-restricted

Despite every directory and file in this project being named "RBP," the
table itself is proteome-scale: 20,431 rows, and confirmed to include
insulin (`P01308`), serum albumin (`P02768`), and hemoglobin alpha
(`P69905`) — none of which are RNA-binding proteins. Known catalogs of
human RBPs list roughly 1,500-2,000 genes; 20,431 is essentially the full
reviewed human proteome. "RBP" describes the lab's research focus, not a
filter on the table's contents. This was asserted incorrectly for a long
stretch of this project's history before being checked directly — worth
remembering when reading anything upstream of this correction that talks
about "the RBP genes" or "RBP-matched" rows.

## Three (at least) competing isoform-table lineages exist in the lab

1. **This project's clean rebuild** — `table_260823_with_rna.csv`,
   verified structurally sound throughout.
2. **The original broken table** —
   `Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv`
   (20.2GB, owned by `fraza`) — the one with duplicate dominant-isoform
   flags and a mostly-empty `UNIQUE` key that motivated this project's
   rebuild in the first place.
3. **`fraza`'s separate April table** —
   `Isoform_Table/Old/Isoform_Post_Merge_PSLab_OpenTargets.csv` (7.36GB,
   different schema again) — what a labmate (`tchhabri`)'s ClinVar variant
   pipeline actually joined against.
4. A labmate (`asharma`) also independently maintains `table_260823.csv`
   and several merged derivatives (`merged_260903A/B.csv`,
   `proteome_with_cdcode_260903.csv`) tracing back to the *same* original
   clean base table this project started from, with condensate-propensity
   predictions (`dG_kT`, `csat_mgmL`, `nu_from_PPII`) and CD-CODE curated
   condensate annotations already merged in.

No single canonical isoform table currently exists across the lab. Before
building anything that depends on "the" isoform table, confirm which
lineage a given file actually traces to.

## The `Dominant_Isoform` flag's actual logic is undocumented anywhere accessible

No script in any directory this project has read access to (`tchhabri`'s
or `fraza`'s) actually *computes* `Dominant_Isoform` — every file either
inherits it wholesale from an upstream table or just reads it. Reverse-
engineered from BRCA1's data: the flagged-dominant row's sequence length
(1,863 aa) matches UniProt's canonical reference length exactly, even
though several *longer* sequences exist among BRCA1's other rows and are
flagged non-dominant. Best guess: it's derived from UniProt's canonical-
isoform designation, not a length- or frequency-based choice — but this
is inference from data, not confirmed from source code. Worth asking
`fraza` directly rather than trusting this write-up as the final answer.

## `tchhabri`'s ClinVar variant-stats table has a confirmed ~10x count inflation

See `07_clinvar_variant_audit/README.md` for the full detail. Short
version: `Total_Pathogenic`/`Total_VUS`/`Total_Benign`, even after properly
deduplicating to one row per protein, sum to 32.7M classified variant
instances — roughly 10x larger than ClinVar's entire actual database.
Likely cause: a missing `VariationID` deduplication step somewhere before
the per-protein rollup (the pipeline's own documentation explicitly warns
this dedup is required). Not yet fixed in that file; this project's
response was to re-pull ClinVar fresh instead (`10_clinvar_repull/`).

## Rare but real CSV corruption in the underlying 315GB ClinVar file

Confirmed via a full 114.9M-row scan (not sampling): a small number of
rows have fields shifted out of alignment — e.g. one row's parsed
`GeneSymbol` field literally contained a fragment of raw JSON text
(`'""review_status"": ""criteria provided'`), and `BRCA1` spuriously
co-occurred with two *other* real genes' ENSGs (`SYNE2`'s and `APC`'s).
Traced to unescaped characters inside the `Germline_Class` JSON column
breaking a naive CSV parser's field boundaries for those specific rows.
Affects only a handful of genes out of tens of thousands — not systemic,
but real. If reusing that file: add a sanity filter (reject rows where
`GeneSymbol` doesn't look like a plausible gene symbol, or where a field
that should be an ID is literally `[]`) rather than trusting every parsed
row at face value.

## Domain-range annotations don't resolve overlaps

See `06_domains_rbd_merge/README.md`. `Classical_RBDs_Annotated.csv`
stores each domain family's range/sequence independently; overlapping
calls from different domain databases are not merged. Confirmed
concretely: one protein's `Domains_total_size` values sum to more
residues than the protein actually has, because one domain call is fully
nested inside another. Fine to use per-family; don't sum across families
without an explicit overlap-merge step first.

## A keyword-name heuristic for biological categorization is directionally right but individually unreliable

See `08_missing_rna_cds_investigation/README.md`. Checking gene/protein
*names* for keywords ("ends in P" → pseudogene, "LINC"/"antisense" →
lncRNA) got the big picture right but mislabeled specific cases — e.g.
tagging a plain pseudogene as "Immunoglobulin gene segment" because its
description happened to mention immunoglobulin-related biology. Always
verify against the authoritative source (Ensembl's own `biotype` field,
in this case) before trusting an individual row's label, even when the
aggregate distribution looks reasonable.

## COSMIC's product names are not self-explanatory about scope

`Cosmic_MutantCensus` *is* restricted to Cancer Gene Census genes (~700-750
genes) as its name suggests, but the **Cancer Mutation Census** file
(`cmc_export.tsv` / `CancerMutationCensus_AllData_*.tsv.gz`) is not — it
covers ~19,996 genes, essentially genome-wide, despite sharing "Census" in
the name. Verified directly (distinct gene count) before trusting either
assumption. Don't infer scope from naming conventions across COSMIC's
product line; check.
