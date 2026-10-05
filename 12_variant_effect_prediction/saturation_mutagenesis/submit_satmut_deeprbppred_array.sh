#!/bin/bash -l
#$ -N deeprbppred_satmut
#$ -cwd
#$ -o /u/home/d/ddcohn/deeprbppred_task_logs/satmut.log.$TASK_ID
#$ -j y
#$ -l h_data=8G,h_rt=1:00:00
#$ -t 1-226
#$ -pe shared 1

module load miniforge/23.11.0
source activate deeprbppred

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_deeprbppred_results
DEEPRBPPRED_DIR=/u/home/d/ddcohn/deeprbppred
WRAPPER=/u/home/d/ddcohn/_claude_deeprbppred_predict_chunk.sh
mkdir -p $OUTDIR

TASK=$(printf "%04d" $SGE_TASK_ID)
IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_TSV=$OUTDIR/chunk_${TASK}.tsv

if [ ! -f "$IN_FASTA" ]; then
  echo "No input chunk for task $SGE_TASK_ID, skipping."
  exit 0
fi

START=$(date +%s)
bash $WRAPPER $IN_FASTA $OUT_TSV $DEEPRBPPRED_DIR balance
RC=$?
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(($(wc -l < $OUT_TSV 2>/dev/null || echo 0) - 1))
echo "Task $SGE_TASK_ID (chunk $TASK): rc=$RC elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$RC" -ne 0 ] || [ "$n_out" -lt "$n_in" ]; then
  exit 1
fi
