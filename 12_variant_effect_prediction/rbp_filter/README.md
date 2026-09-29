# Filtering the variant-effect tables down to RNA-binding proteins

The SpliceAI/DeepLoc/protGPS pipeline in `../splicing/` and
`../localization_condensate/` was run genome-wide -- on every GRCh38
ClinVar variant and every usable CMC variant, regardless of gene. That's
correct for the pipeline itself (nothing about SpliceAI, DeepLoc, or
protGPS is RBP-specific), but this project's actual scope is RNA-binding
proteins specifically, so the six merged tables need to be filtered down
before use, not read as-is.

## RBP gene list

The active filter list, `rbp_gene_symbols.txt` (2,215 genes), is the
**union of two independent RBP definitions**:

1. **Domain-based (533 genes, `rbp_gene_symbols_domain_based.txt`).**
   `Has_RBD == 1` in `/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv`
   (this project's own gene-level isoform table) -- proteins with an
   InterPro-annotated classic RNA-binding domain (RRM, KH, DRBM,
   helicase, SAM, PAZ/Piwi, La-motif, etc; 33 domain names total). That
   table has no gene-symbol column, so `build_rbp_gene_list.py` maps its
   534 NCBI gene IDs to symbols via the raw ClinVar file's own
   `GeneID`/`GeneSymbol` columns (533/534 resolved; one gene ID,
   `100913187`, has no ClinVar match and is a known, negligible gap).
   This definition is narrow by construction -- it misses RBPs that bind
   RNA via disordered regions, non-canonical folds, or domains outside
   that 33-name list (e.g. zinc fingers, LSm/Sm, cold-shock, RGG motifs
   are notably absent from it).

2. **RBP2GO high-confidence human RBPs (1,968 genes,
   `rbp_gene_symbols_rbp2go.txt`).** Downloaded from
   [RBP2GO](https://rbp2gov2.dkfz.de/) (References & Data -> Protein
   Data -> Homo sapiens -> RBP Dataset;
   `High_confidence_human_RBPs_rbp2go.txt`, 1,977 UniProt IDs with
   `RBP_status == RBP`), an aggregation of many RNA-interactome-capture
   (CLIP-seq/mass-spec) studies, so it also captures RBPs with no
   annotated domain at all. `resolve_rbp2go_uniprot_ids.py` maps those
   UniProt IDs to gene symbols via the labmate's isoform table
   (`gene_symbol` + `uniprot_secondary_accessions` +
   `swissprot_canonical_accessions` columns), resolving 1,968/1,977 (9
   obscure/fragment UniProt IDs didn't match anything).

**Overlap between the two: only 286 genes.** 247 domain-based genes
aren't in RBP2GO's high-confidence set (plausible domain-present-but-
unconfirmed cases, not investigated further), and 1,682 RBP2GO genes
have no classic RBD (exactly the unconventional-RBP gap the domain-only
list was expected to miss). `build_union_list.py` takes the union of
both (2,215 genes) rather than picking one -- the decision made was to
keep any gene flagged by either method, on the reasoning that the
247 domain-based-only genes might be real RBPs RBP2GO's aggregated
studies simply haven't captured yet, rather than risk dropping them.

Output: `rbp_gene_symbols.txt` (2,215, the active list used by
`filter_rbp_tables.py`), plus the two source lists it was built from,
all checked into this repo since they're small and are the actual
filter definitions used -- anyone rerunning `filter_rbp_tables.py`
reproduces the same filter directly, without recomputing the ID
resolution steps.

## Filtering the six tables

`filter_rbp_tables.py` filters each of the six merged tables in
`/u/project/kappel/ddcohn/protein_variant_effects/` down to rows whose
gene is in `rbp_gene_symbols.txt`, writing a `*_rbp.tsv` sibling next to
each original (not checked into git -- covered by the same data
`.gitignore` rule as the unfiltered tables). DeepLoc/protGPS deltas
already carry a gene column directly; ClinVar's SpliceAI table and CMC's
SpliceAI table don't (SpliceAI needs no gene/protein info to run), so
those are filtered via a join: ClinVar through
`clinvar_metadata_lookup.tsv` (`VariationID` -> `GeneSymbol`, built from
every GRCh38 row), and CMC through a fresh `GENOMIC_MUTATION_ID` ->
`GENE_NAME` map read directly from the raw CMC file (`cmc_protein_changes.tsv`
only covers the protein-editable variant categories, not the
synonymous/frameshift/noncoding variants SpliceAI also scored, so it
isn't a complete enough join key on its own).

## Row counts, before -> after

| Table | Genome-wide | RBP-only (2,215 genes) |
|---|---|---|
| ClinVar DeepLoc deltas | 2,657,205 | 259,472 |
| ClinVar protGPS deltas | 2,587,382 | 253,127 |
| ClinVar SpliceAI scores | 4,137,146 | 388,515 |
| ClinVar metadata lookup | 4,494,430 | 426,887 |
| CMC DeepLoc deltas | 4,283,459 | 381,249 |
| CMC protGPS deltas | 4,190,344 | 372,819 |
| CMC SpliceAI scores | 5,070,804 | 450,667 |

(An earlier pass filtered on the 533-gene domain-based list alone; those
numbers are superseded by the union-list numbers above. See git history
for the domain-only row counts if needed for comparison.)

## What's not done yet

The notebooks in `../../notebooks/` still read the unfiltered,
genome-wide tables. RBP-only versions of the exploration notebooks
haven't been built yet.
