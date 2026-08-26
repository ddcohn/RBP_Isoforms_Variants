# Isoform Table diagnostics (old, broken table only)

Ran against the original 316-column, 68,645-row table
(`Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv`) to
characterize why it looked wrong. Not part of the rebuild pipeline — kept for
reference on what was broken and why the table was rebuilt from a clean base
instead of patched.

- **`sanity_check.py`** — first-pass structural check: row counts, `row_kind`
  distribution, `dominant_isoform` value distribution, duplicate `UNIQUE`
  values, uniprot_id/protein_key cardinality, ENSG/ENSP emptiness. Found the
  core problem: `dominant_isoform` is just a recoded copy of `row_kind`
  (`swissprot_canonical` → 1), not a real one-per-protein flag — 3,872
  proteins had more than one row flagged dominant.

- **`no_dominant_check.py`** — follow-up check specifically for proteins with
  *zero* rows flagged dominant (534 found, including some with dozens of
  isoform rows and none flagged, e.g. `ABI2` with 88 rows).

- **`completeness_check.py`** — cross-references a checklist of expected
  columns (from the original project spec) against what's actually populated
  in the table, to find columns that are present but mostly empty (e.g.
  condensate data only 9.3% populated) versus genuinely absent (e.g.
  `canonical RBDs`, granular OpenTargets evidence fields).
