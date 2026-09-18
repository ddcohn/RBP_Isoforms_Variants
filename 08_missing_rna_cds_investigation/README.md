# Why 876 proteins are missing RNA/CDS sequence — categorized, then verified

`../02_rna_cds_pipeline/exact_match_rebuild.py` leaves 876 proteins (4.3%)
without a resolved RNA/CDS sequence. This directory categorizes *why*,
then — importantly — checks whether that categorization is actually
correct rather than trusting a quick heuristic.

## Scripts, in run order

1. **`categorize_missing_rna.py`** — for each of the 876 proteins, pulls
   the UniProt protein name and gene name(s), then sorts into biological
   categories using **keyword matching** on those names (e.g. "pseudogene"
   if the gene symbol ends in P, "lncRNA" if the name contains "antisense"
   or "LINC", etc.). Fast, but a heuristic — see the verification step
   below for how well it actually held up.

2. **`plot_initial_distribution.py`** — bar chart of the keyword-heuristic
   category counts. Superseded by the corrected version below; kept to
   show what changed and by how much.

3. **`verify_categories_against_ensembl.py`** — re-checks every one of the
   598 rows that have an `ensembl_gene` against Ensembl's **authoritative**
   gene biotype (not name guessing). Findings:
   - The core conclusion holds: **zero** of the 876 turned out to be
     genuinely protein-coding genes wrongly excluded.
   - But the keyword heuristic's breakdown was unreliable: only 57.5% of
     rows agreed with the real biotype once normalized. Most disagreement
     (172 + 73 of 254 mismatches) was the vague "Other/unclassified"
     catch-all getting correctly resolved into "Pseudogene" or
     "Antisense/lncRNA" — not wrong, just imprecise.
   - **12 rows had a genuinely wrong, specific label** — e.g. `FCGR1BP`
     tagged "Immunoglobulin gene segment" and `TAAR3P` tagged "TCR gene
     segment" by loose substring matching, when both are actually plain
     pseudogenes unrelated to Ig/TCR.
   - The "No ensembl_gene at all" bucket was itself undercounted: only 140
     of the true 278 gene-ID-less rows were labeled that way; the other
     138 got swept into name-based guesses despite having no gene ID to
     search from in the first place — a script-ordering bug (name
     heuristic ran before checking whether an ID even existed).

4. **`plot_corrected_distribution.py`** — bar chart of the real, verified
   distribution. Notably, **Pseudogene becomes the largest category
   (35.5%, not 16.2%)**, and the Ig/TCR counts drop from 16/6 to 2/1 once
   corrected.

## Corrected distribution (of 876)

| Category | Count | % |
|---|---|---|
| Pseudogene | 311 | 35.5% |
| Antisense/lncRNA/uncharacterized locus | 267 | 30.5% |
| No `ensembl_gene` at all | 278 | 31.7% |
| Ensembl ID not found / retired | 9 | 1.0% |
| Mitochondrial-encoded micropeptide | 8 | 0.9% |
| Immunoglobulin gene segment | 2 | 0.2% |
| TCR gene segment | 1 | 0.1% |

**Takeaway for anyone reusing this pattern**: a fast name-keyword heuristic
is fine for a rough first pass, but individual category labels from it
should not be trusted without checking against the authoritative source
(here, Ensembl's own `biotype` field) — especially for any category that
implies something narrow and specific (Ig/TCR/mitochondrial), where a
loose substring match is easy to get wrong.
