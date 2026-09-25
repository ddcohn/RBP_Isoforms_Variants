# Subcellular localization and condensate propensity (DeepLoc, protGPS)

Two predictors, run on both the wild-type and mutant sequences from
`../mutant_sequence_construction/` (2,657,205 mutant + 19,405 WT), so the
effect of a variant is the *change* in predicted score, not the
wild-type score alone:

- **DeepLoc 2.1** — general subcellular localization (nucleus, cytoplasm,
  membrane, etc.), via ESM1b (650M params).
- **protGPS** — biomolecular condensate localization specifically
  (stress granule, nucleolus, P-body, etc. — a narrower, condensate-
  focused prediction than DeepLoc's general compartments), via ESM2 (8M
  params — the small ESM2 checkpoint, not the large one, and no DR-BERT
  needed despite that being unclear from protGPS's own repo structure
  going in).

Adapting protGPS required reading its `notebook/Predict.ipynb` directly
(no clean CLI) — `protgps_predict_validation.py` is that logic pulled out
into a standalone script, validated against the notebook's own documented
example output (exact match) before running at scale.

## Real throughput numbers (not estimates) — why this needed the shared GPU pool

A 1000-sequence timing test on the Kappel lab's dedicated GPU node
(2x L40S) gave:

- **DeepLoc: ~4.4 seqs/sec/GPU.** At that rate, the full ~2.68M sequences
  would take **~3.5 days even split across both dedicated GPUs** — by far
  the bottleneck of the two tools.
- **protGPS: ~49 seqs/sec/GPU** at batch_size=8 — full set in ~7.5hr
  across both GPUs. batch_size=32 crashed with a CUDA OOM, not from GPU
  contention (confirmed `CUDA_VISIBLE_DEVICES` was correctly pointing at
  the idle GPU) but from real memory pressure: sequence lengths in this
  dataset have a long tail (median 857 residues, but max 36,003 —
  titin-scale), and attention memory scales with the *longest* sequence
  in a batch, not the average.

To get DeepLoc down from days to hours, both tools were instead run as
20-task array jobs (`build_deeploc_chunks.py`, reused for both since they
take the same input) spread across Hoffman2's shared GPU pool
(`gpu_a100.q`/`gpu_l40s.q`/`gpu_h100.q`/`gpu_RTX2080Ti.q`, plus others
only visible once actually requested: `gpu_smp.q`/`gpu_rh7.q`/
`gpu_v100.q`), not just the 2 dedicated GPUs — see the top-level
`12_variant_effect_prediction/README.md` for the resource-flag details
that make a job eligible for the shared pool vs. the dedicated queue.

`protgps_predict_scale.py` additionally sorts sequences by length before
batching (minimizes padding waste and stops one long outlier from
blowing up an otherwise-short batch), adapts batch size down for longer
sequences (2 for >2000 residues, 4 for >1000, 8 otherwise), skips
sequences over 5,000 residues entirely (~2.5% of the set, flagged `NA`
in the output rather than silently dropped) as a safety margin against
unpredictable OOM crashes on the longest outliers, and falls back to
one-at-a-time processing with a per-sequence OOM guard so one bad
sequence can't take down an entire chunk's worth of otherwise-successful
predictions.

## Status

| Set | Tool | Result |
|---|---|---|
| ClinVar (20 chunks) | DeepLoc 2.1, protGPS | complete |
| COSMIC CMC (32 chunks) | protGPS | complete |
| COSMIC CMC (32 chunks) | DeepLoc 2.1 | complete (chunks 009, 017 rerun via `submit_deeploc_cmc_retry.sh`; 026 via a backup run) |

CMC scripts: `build_cmc_deeploc_chunks.py`, `submit_cmc_deeploc_array.sh`,
`submit_cmc_protgps_array.sh`, `protgps_predict_cmc_scale.py`,
`merge_cmc_protgps_deltas.py`. The CMC mutant sequences come from
`../mutant_sequence_construction/` (`parse_cmc_protein_changes.py`,
`fetch_cmc_wt_sequences_grch37.py`, `build_cmc_mutant_sequences.py`).
DeepLoc-to-delta merging and the DeepLoc notebooks are not written yet.

## Lessons from the runs (read before rerunning)

- **Exit codes:** the original array scripts ended with an `echo`, so a
  task that crashed (CUDA illegal memory access on CMC chunk 017) was
  recorded by SGE as exit 0. The scripts now capture `RC=$?` and
  `exit $RC`. Even so, `qacct`/`qstat` are not evidence of success.
- **Verify by content:** a DeepLoc results CSV must have
  (sequences in the chunk FASTA + 1 header) lines; protGPS TSVs likewise.
  Checking that files exist or are non-empty missed truncated outputs.
- **Duplicate IDs:** CMC chunk 017 contains one FASTA ID
  (`COSV52782716_ENST00000269305`) twice with two different sequences.
  DeepLoc keeps one row per ID, so that chunk's CSV has 134,484 lines
  for 134,484 sequences (one fewer than the rule above): one variant has
  no DeepLoc prediction. It is the only duplicated ID in the CMC set.
- **GPU choice:** DeepLoc takes ~2 h per 134k-sequence chunk on an L40S
  but 20 h or more on the P4 GPU (node m3), so pin to L40S/A100/H100
  queues. Shared GPU queues cap `h_rt` at 24 h and need
  `-l cuda=1,require_gpu=1` (no `highp`); `kappel_gpu.q` needs `highp`.
  The output directory's parent must exist before `deeploc2` runs.
