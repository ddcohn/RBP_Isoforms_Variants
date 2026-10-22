#!/bin/bash
# Syncs this repo's pipeline code (not data) to the canonical Hoffman2
# location the Snakefile expects. Run this locally after any script
# change, before running the pipeline on the cluster -- the Snakefile's
# paths all resolve relative to where IT lives, so the whole tree below
# has to move together, not just the one script you edited.
set -euo pipefail

LOCAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_HOST="ddcohn@hoffman2.idre.ucla.edu"
REMOTE_DIR="/u/project/kappel/ddcohn/RBP_pipeline_scripts"

ssh "$REMOTE_HOST" "mkdir -p $REMOTE_DIR"

for dir in workflow saturation_mutagenesis rbp_filter localization_condensate \
           mutant_sequence_construction splicing; do
  rsync -avz --exclude ".DS_Store" --exclude "__pycache__" \
    "$LOCAL_DIR/$dir/" "$REMOTE_HOST:$REMOTE_DIR/$dir/"
done

echo "Deployed to $REMOTE_HOST:$REMOTE_DIR"
echo "Run the pipeline with:"
echo "  ssh $REMOTE_HOST"
echo "  module load miniforge/23.11.0 && source activate snakemake"
echo "  cd $REMOTE_DIR/workflow"
echo "  snakemake --profile profiles/hoffman2 --cores 1 all"
