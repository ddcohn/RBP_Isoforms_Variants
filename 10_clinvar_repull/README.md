# Re-pulling ClinVar from scratch

The lab's existing ClinVar-derived variant table (audited in
`../07_clinvar_variant_audit/`) had two confirmed problems: pathogenic/VUS/
benign counts inflated ~10x, and rare but real CSV corruption. Rather than
untangle which step in that file's long processing chain (62GB XML release
→ isoform join → repeated `fix_`/`verify_`/`reannotate` patches) caused
which issue, this re-pulls ClinVar directly from NCBI's own bulk file —
same "clean rebuild over patching" approach as the original isoform table.

## What we checked before pulling anything

NCBI publishes `variant_summary.txt` — a flat, official, tab-delimited
report, one row per variant. Before switching to it, we verified against
NCBI's own documentation (not assumed):
- **No variants lost** — same full ClinVar variant universe as the XML
  release used before, not a subset.
- **A different duplication axis to handle**: one row per variant **per
  genome assembly** — most variants appear twice (GRCh37 and GRCh38).
  Confirmed directly (`AlleleID` 15041 appears twice, identical except
  `Assembly`/coordinates). Filter to one assembly before counting anything.
- **A real trade-off**: this file *aggregates* multiple submitters'
  classifications into summary fields (`ClinicalSignificance`,
  `ClinSigSimple`, `NumberSubmitters`) rather than keeping each submitter's
  individual call, condition, and review status separate (which the old
  file did, as a per-variant JSON array). Traded submitter-level provenance
  for a simpler, more reliable structure.

## `extract_clinvar_grch38.py`

Downloads `variant_summary.txt.gz` from NCBI's public FTP (443MB
compressed, no credentials needed), filters to `Assembly == GRCh38`, and
writes out a focused column set: `VariationID` (dedup key), `GeneSymbol`/
`GeneID`/`HGNC_ID`, `ClinicalSignificance`/`ClinSigSimple`, `ReviewStatus`/
`NumberSubmitters`, `PhenotypeList`, `Origin`/`OriginSimple`, position/
allele fields, and `SomaticClinicalImpact`/`Oncogenicity`.

Result: 4,494,430 GRCh38 rows (confirmed against a separate `awk` count of
the raw file before running the Python extraction — matched exactly).
4,492,073 distinct `VariationID`s (99.95% unique) — the small residual
(2,357 duplicates) is almost certainly legitimate cases where one ClinVar
"Variation" spans multiple `AlleleID` records, not a data problem. This
scale (~4.5M distinct variants) matches ClinVar's actual known database
size, unlike the ~10x-inflated numbers from the old derived file — a good
sanity check that this fresh pull isn't repeating the same problem.

**Note**: `variant_summary.txt` has no `ENSEMBL_ID` column — only
`GeneID` (Entrez) and `HGNC_ID`. Joining this to the isoform table (which
has `ensembl_gene` directly) needs an extra hop through a mapping file —
`Cosmic_Genes` (see `../09_cosmic_pipeline/`) already has `ENTREZ_ID`,
`HGNC_ID`, and `GENE_ACCESSION` (ENSG) together in one row and can serve
as that bridge, avoiding a second lookup.

## Not yet built

- The actual join to the isoform table / to the fresh COSMIC pull (matching
  by gene, per the decision above — ClinVar and COSMIC agreed to land as
  two separate tables, joined by gene when needed, not merged into one).
- Whatever downstream annotation (SASA, domain/disorder overlap) the old
  file had layered on top — not recomputed here.
