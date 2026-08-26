# RNA / CDS sequence pipeline

Adds `rna_sequence` (full spliced mRNA) and `cds_sequence` (coding region
only) columns to the clean base table, keyed by `uniprot_accession`.

## Evolution of the approach (read in this order to understand *why*)

1. **`pull_rna.py`** — first attempt: map `uniprot_accession` →
   `Ensembl_Transcript` via UniProt's bulk ID-mapping API, take the first
   transcript returned per protein. Works, but "first returned" isn't
   necessarily the transcript matching the specific isoform on file — many
   proteins map to multiple transcripts (multi-isoform genes), so this
   silently picks the wrong protein a meaningful fraction of the time.

2. **`pull_rna_v2.py`** — the real fix: use UniProt's `uniprotkb/accessions`
   bulk endpoint with `fields=accession,xref_ensembl`, which returns each
   Ensembl transcript cross-reference tagged with the *specific UniProt
   isoform* it encodes (`isoformId`, e.g. `P04637-1`). Only accept a
   transcript whose isoform tag matches the canonical accession (no suffix,
   or explicit `-1`). Falls back to "first available" only when no
   canonical-tagged transcript exists (`transcript_selection_method` column
   records which case applied: `canonical_isoform_match`,
   `fallback_no_canonical_tag`, `no_ensembl_xref`).
   **Note:** Ensembl's `/sequence/id` endpoint rejects versioned transcript
   IDs (`ENST00000424496.3`) — strip the version before fetching sequence.

3. **`patch_missing_rna.py`** — retries transcripts that got a row in
   `pull_rna_v2.py`'s output but no sequence (transient Ensembl 5xx errors
   that weren't retried the first time). Adds retry-on-503.

4. **`final_patch.py`** — handles the last few unrecoverable cases: some
   UniProt cross-references point to Ensembl transcript IDs that have been
   permanently retired (confirmed via `/archive/id/:id`, which shows
   `possible_replacement: []`). For these, manually tries alternate
   transcripts tagged to the same isoform; a couple of proteins are
   confirmed truly unrecoverable this way.

5. **`categorize_no_map.py`** — categorizes the ~1,030 proteins with *no*
   Ensembl protein-coding transcript cross-reference at all. They fall into
   explainable biological categories: pseudogenes (frameshifted/disrupted,
   e.g. `CYP2D7`), immunoglobulin/TCR gene segments (assembled by somatic
   DNA recombination, no fixed reference sequence), and micropeptides
   translated from a short internal ORF inside a gene Ensembl otherwise
   classifies as non-coding RNA (e.g. `PINT87aa` inside the lncRNA
   `LINC-PINT`, or humanin-like peptides inside the mitochondrial 16S rRNA
   gene `MT-RNR2`). None of these are pipeline bugs.

6. **`pull_cds.py`** — separate pass to also fetch `cds_sequence` (type=cds)
   for every transcript already selected as `ensembl_transcript_used`,
   reusing the same transcript choice as the RNA pull (no re-resolution
   needed).

7. **`categorize_fixable.py`** — re-examines the ~1,030 no-RNA proteins,
   this time checking whether the *base table's own* `ensembl_gene` column
   (populated independently of the UniProt cross-reference used above) has a
   value. 752 do. Categorizes those into the same biological buckets as
   above, isolating a genuinely "unclassified, might be fixable" subset
   (362) — later shown to still be mostly the same categories (pseudogenes,
   retrogenes) just missed by naming-convention heuristics.

8. **`gene_level_fallback.py`** — attempts to recover RNA for the 752 rows
   above by looking up the gene directly (`/lookup/id/{ENSG}?expand=1`) and
   accepting *any* transcript with a `Translation` object (not filtering by
   `biotype == protein_coding`, since Ensembl tags real translated
   transcripts for e.g. immunoglobulin constant genes as `IG_C_gene`, not
   `protein_coding`). Recovered 156 rows, **but only ~9% were verified
   correct** (see `03_sanity_checks/`) — this method has no isoform
   verification, so it frequently grabs the wrong transcript for genes with
   multiple protein products. Superseded by `exact_match_rebuild.py`.

9. **`exact_match_rebuild.py`** — the principled fix for both the
   `gene_level_fallback` weakness and an earlier discovery that
   `fallback_no_canonical_tag` rows (401 of them) had the exact same
   problem. For every candidate protein, fetches **every** transcript of its
   gene (`ensembl_gene`), translates each one's CDS, and keeps only the
   transcript whose translation is a **byte-for-byte exact match** to the
   stored protein sequence — no guessing, no "canonical" assumption.
   Verified 25/25 exact matches on a test sample. Has JSON checkpointing at
   every phase (gene lookup, CDS fetch, mRNA fetch) since this is a
   long-running job (~8+ hours for gene lookups alone at observed API
   throughput of ~0.3–1.9s/gene, which varies with Ensembl server load).
   **Must be run via `qsub`, not directly on the login node** — see below.

10. **`submit_exact_match.sh`** — UGE/SGE job script for `exact_match_rebuild.py`.
    Requests `h_data=8G,h_rt=24:00:00`. Needed because the Hoffman2 login
    node enforces a ~1GB memory / 1hr CPU `ulimit` on background processes
    and silently kills anything that exceeds it (no error, no traceback —
    the process just vanishes). Submit with:
    ```
    /u/local/bin/qsub submit_exact_match.sh
    /u/local/bin/qstat -u <username>   # check status
    ```
    If the job hits its `h_rt` wall-clock limit before finishing, just
    resubmit the same script — checkpoints let it resume rather than restart.

## Column reference

| Column | Meaning |
|---|---|
| `rna_sequence` | Full spliced mRNA (5'UTR + CDS + 3'UTR) |
| `cds_sequence` | Coding region only |
| `ensembl_transcript_used` | Which `ENST` the sequences came from |
| `transcript_selection_method` | How it was chosen — see script list above for what each value means |
