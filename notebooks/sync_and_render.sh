#!/bin/bash
# Pulls the 7 tables the Hoffman2 pipeline produces
# (12_variant_effect_prediction/workflow/, target `all`) down to this
# directory and re-executes all 6 exploration notebooks against them.
# Run from notebooks/. Requires an ssh alias/host for Hoffman2 as $1
# (e.g. ddcohn@hoffman2.idre.ucla.edu).
set -euo pipefail

HOST="${1:?usage: sync_and_render.sh <ssh-host>}"
REMOTE=/u/project/kappel/ddcohn/protein_variant_effects

FILES=(
  clinvar_spliceai_scores.tsv cmc_spliceai_scores.tsv
  clinvar_deeploc_deltas.tsv cmc_deeploc_deltas.tsv
  clinvar_protgps_deltas.tsv cmc_protgps_deltas.tsv
  clinvar_metadata_lookup.tsv
)

for f in "${FILES[@]}"; do
  scp -q "${HOST}:${REMOTE}/${f}" .
done

python3 run_notebook.py \
  spliceai_clinvar_exploration.ipynb spliceai_cmc_exploration.ipynb \
  deeploc_clinvar_exploration.ipynb deeploc_cmc_exploration.ipynb \
  protgps_clinvar_exploration.ipynb protgps_cmc_exploration.ipynb
