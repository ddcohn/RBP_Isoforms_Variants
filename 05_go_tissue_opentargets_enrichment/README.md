# GO terms, tissue expression, and OpenTargets disease associations

Adds three more enrichment layers to the clean base table, keyed by
`uniprot_accession` / `ensembl_gene`: Gene Ontology terms, bulk tissue
expression, and disease associations.

## Scripts, in run order

1. **`pull_go.py`** — pulls Cellular Component, Biological Process, and
   Molecular Function GO terms per protein from UniProt's REST API
   (`fields=go_p,go_c,go_f`), batched 100 accessions per request. Stores each
   aspect as a semicolon-separated list of `GO_ID:term_name` pairs.
   Coverage: 92.8% CC, 85.4% BP, 80.0% MF.

2. **`go_result_to_csv.py`** — flattens the raw per-protein GO JSON result
   into a standalone CSV for quick inspection before merging.

3. **`pull_tissue.py`** — pulls **bulk** tissue expression from OpenTargets'
   `baselineExpression` field, explicitly filtering out per-cell-type
   single-cell breakdowns (Tabula Sapiens) by keeping only rows where
   `celltypeBiosampleFromSource` is empty. Without this filter the full
   single-cell breakdown is ~29M rows across the whole gene set — the
   filtered bulk version is ~2.8M rows, comparable in scale to the other
   pulls. **Note:** OpenTargets' GraphQL API rejects a query as "too
   expensive" above roughly `batch_size × page_size ≈ 10,000` — this script
   uses 5 genes/request × 1,500 rows/gene (empirically, even
   heavily-studied genes cap out around ~1,450 total baseline-expression
   rows, so 1,500 is safe headroom).

4. **`merge_all_pending.py`** — joins the GO, tissue, and OpenTargets disease
   association results (the latter computed earlier in the original
   pipeline but deliberately *not* merged in at the time — see comment in
   `pull_opentargets.py` in `../02_rna_cds_pipeline/`, which held off merging
   because `exact_match_rebuild.py` might still have been running against
   the same target file) onto the base table in one pass.

5. **`verify_merge.py`** — the integrity check run after *every* merge in
   this project from this point on: confirms row count is unchanged, no
   ragged rows, every pre-existing column is byte-identical between the
   pre-merge backup and the merged file, and reports coverage of the newly
   added columns. Always back up the target file before merging and run
   this after.

`submit_go.sh` / `submit_tissue.sh` / `submit_merge.sh` / `submit_verify.sh`
are the corresponding UGE job scripts — the GO and tissue pulls are
long-running API-batch jobs, and the merge/verify steps need enough memory
to hold the (by then hundreds-of-MB) table in memory twice.

## Column reference

| Column | Meaning |
|---|---|
| `GO_Cellular_Component` / `_Biological_Process` / `_Molecular_Function` | `;`-separated `GO_ID:term` pairs |
| `Tissue_Expression_Bulk` | `;`-separated `tissue\|datasourceId\|median\|unit` entries |
| `OpenTargets_DiseaseId` / `_DatatypeId` / `_Score` | parallel `;`-separated lists, one entry per disease association |
