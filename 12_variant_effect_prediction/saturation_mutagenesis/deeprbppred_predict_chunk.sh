#!/bin/bash
# Run Deep-RBPPred on one chunk FASTA and pair its per-line scores back with
# the mutation_id headers (Deep-RBPPred.py only prints raw scores, one per
# line, in input order -- it doesn't echo back sequence identity itself).
set -euo pipefail

IN_FASTA=$1
OUT_TSV=$2
DEEPRBPPRED_DIR=$3   # directory containing Deep-RBPPred.py, feature, RBP_Model/
MODEL=${4:-balance}

cd "$DEEPRBPPRED_DIR"
grep '^>' "$IN_FASTA" | sed 's/^>//' > /tmp/deeprbppred_ids_$$.txt
python Deep-RBPPred.py "$IN_FASTA" "$MODEL" > /tmp/deeprbppred_scores_$$.txt 2>/tmp/deeprbppred_stderr_$$.txt

n_ids=$(wc -l < /tmp/deeprbppred_ids_$$.txt)
n_scores=$(wc -l < /tmp/deeprbppred_scores_$$.txt)
if [ "$n_ids" != "$n_scores" ]; then
  echo "MISMATCH: $n_ids ids vs $n_scores scores -- not writing output" >&2
  cat /tmp/deeprbppred_stderr_$$.txt >&2
  rm -f /tmp/deeprbppred_ids_$$.txt /tmp/deeprbppred_scores_$$.txt /tmp/deeprbppred_stderr_$$.txt
  exit 1
fi

echo -e "mutation_id\trbp_probability" > "$OUT_TSV"
paste /tmp/deeprbppred_ids_$$.txt /tmp/deeprbppred_scores_$$.txt >> "$OUT_TSV"

rm -f /tmp/deeprbppred_ids_$$.txt /tmp/deeprbppred_scores_$$.txt /tmp/deeprbppred_stderr_$$.txt
