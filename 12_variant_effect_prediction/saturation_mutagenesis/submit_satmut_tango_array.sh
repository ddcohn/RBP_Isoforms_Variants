#!/bin/bash -l
#$ -N tango_satmut
#$ -cwd
#$ -o /u/home/d/ddcohn/tango_task_logs/satmut.log.$TASK_ID
#$ -j y
#$ -l h_data=2G,h_rt=8:00:00
#$ -t 1-226
#$ -pe shared 1

module load miniforge/23.11.0

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_tango_results
TANGO_BIN=/u/home/d/ddcohn/tango_x86_64_release
DRIVER=/u/home/d/ddcohn/_claude_tango_predict_chunk.py
mkdir -p $OUTDIR

TASK=$(printf "%04d" $SGE_TASK_ID)
IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_TSV=$OUTDIR/chunk_${TASK}.tsv

if [ ! -f "$IN_FASTA" ]; then
  echo "No input chunk for task $SGE_TASK_ID, skipping."
  exit 0
fi

START=$(date +%s)
python3 $DRIVER $IN_FASTA $OUT_TSV $TANGO_BIN
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(($(wc -l < $OUT_TSV 2>/dev/null || echo 0) - 1))
echo "Task $SGE_TASK_ID (chunk $TASK): elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$n_out" -lt "$n_in" ]; then
  echo "Incomplete -- resubmit this task to resume from where it left off."
  exit 1
fi
