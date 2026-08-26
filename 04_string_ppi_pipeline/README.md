# STRING protein-protein interaction pipeline

Adds `PPI_UniProt_Partners` (all interaction partners, as UniProt accessions)
and `PPI_UniProt_Partners_in_Dataframe` (the subset also present in this
dataset), computed three independent ways for cross-validation.

## Two STRING quirks that will silently corrupt results if missed

1. **Default result cap.** STRING's `interaction_partners` endpoint silently
   returns only the top 10 partners per protein unless you pass an explicit
   `limit` override (e.g. `limit=100000`). A well-studied protein like FUS
   actually has ~90-160 high-confidence partners depending on STRING
   version — the un-overridden default would make it look like it has 10.

2. **Silent ID renaming.** STRING normalizes/resolves some input protein IDs
   to a different current canonical ID before returning interaction data.
   Querying with the original ID and expecting the response's `stringId_A`
   to match it will silently drop results for any protein STRING renamed.
   Fix: call `get_string_ids` first to resolve every input ID, then query
   `interaction_partners` using the *resolved* ID and match results back
   against that.

## Scripts

- **`pull_ppi.py`** — the full pipeline against STRING's live REST API
  (`string-db.org`, confirmed to be serving v12.0 — check
  `version.string-db.org/api/json/version`, don't assume it's current):
  resolve IDs → query partners (`required_score=700`, "high confidence",
  no cap) → map partner ENSPs to UniProt accessions via UniProt's bulk ID
  mapping → write both columns. Has checkpointing and auto-retry/resubmit
  for the UniProt mapping job (which had a real transient backend outage
  during one run — "Connection refused" from UniProt's own infrastructure,
  not a bug in this code).

- **`compare_cds.py`** — separate investigation: for the ~13,378 proteins
  that mapped to multiple Ensembl transcripts, do all their transcripts
  encode the same coding sequence, or different ones? (Answer needed before
  trusting "first transcript" as a selection strategy elsewhere in the
  pipeline.) Long-running (92,341 transcripts to fetch); has checkpointing.

### STRING v12.5 comparison (the live API turned out to be on v12.0)

- **`resolve_v125.py`** — resolves all query protein IDs against STRING
  v12.5's version-pinned endpoint (`version-12-5.string-db.org` — its own
  `/version` field misreports "12.0", but the actual interaction data is
  confirmed v12.5 by cross-checking against a locally downloaded STRING
  v12.5 flat file). Saved once and shared between the two v12.5 methods
  below so they're a true apples-to-apples comparison.

- **`count_partners.awk`** / **`count_partners_resolved.awk`** — single-pass
  AWK scripts that scan a locally downloaded STRING `protein.links` flat
  file (`9606.protein.links.v12.5.txt`, ~10.8M rows) and extract
  high-confidence pairs for a given set of query proteins. Far faster than
  API calls (whole-file scan in ~11s) once the file is downloaded locally.
  The `_resolved` variant takes STRING-resolved IDs as input (see above) —
  the non-resolved version undercounts because it can't match STRING's
  internally-renamed IDs.

- **`finish_local_v125.py`** / **`finish_local_v125_resolved.py`** —
  complete the local-file method by mapping partner ENSPs to UniProt
  accessions. **Note the bug fixed between v1 and v2**: partner ENSPs from
  the flat file carry a `9606.` species prefix, but UniProt's
  `Ensembl_Protein` ID mapping expects bare IDs — submitting prefixed IDs
  silently returns zero matches (job "succeeds" but maps nothing).

- **`pull_ppi_v125_api.py`** — the live-API equivalent, using the same
  pre-resolved IDs from `resolve_v125.py` for direct comparability.

- **`merge_v125_columns.py`** — merges both v12.5 methods into the table as
  four new columns (`PPI_UniProt_Partners_v125_local`,
  `PPI_UniProt_Partners_in_Dataframe_v125_local`, and `_api` equivalents),
  alongside the original v12.0 columns — all three kept side by side rather
  than overwritten, so they can be compared directly.

## Validation

- **`sanity_check_ppi.py`** — structural checks (self-references, whether
  the in-dataframe list is a subset of the full partner list, duplicates,
  malformed accessions), partner-count distribution, and a symmetry check
  (if A lists B as a partner, does B list A back?). Found 78 legitimate
  self-references, all traced to multi-copy gene families (e.g. histones)
  where UniProt collapses several paralogous genes into one accession.
  93.9% of in-dataframe edges are reciprocated; the rest are attributable to
  proteins with multiple Ensembl IDs where only the first was queried.

- **`compare_methods_per_row.py`** — row-by-row agreement, not just
  aggregate counts: Jaccard similarity between the two v12.5 methods
  (same STRING version, different access method — validates the pipeline
  itself: 49.8% exact match, 75.9% Jaccard ≥0.9) versus v12.0-vs-v12.5-API
  (different STRING version — measures real database growth between
  releases, expected to differ more: median Jaccard 0.600).

## Version note

STRING's live default API (`string-db.org`) was confirmed serving v12.0 at
the time this was built, while a locally downloaded flat file and the
version-pinned `version-12-5.string-db.org` endpoint have the newer v12.5
data (~69% more partners for some proteins, e.g. FUS: 92 → 160-162). Check
`https://version.string-db.org/api/json/version` before assuming which
version you're querying.
