#!/bin/bash -l
#$ -N protgps_satmut
#$ -cwd
#$ -o /u/home/d/ddcohn/protgps_task_logs/satmut.log.$TASK_ID
#$ -j y
#$ -q kappel_gpu.q
#$ -l h_data=24G,h_rt=2:00:00,highp,require_gpu=1,cuda=1
#$ -t 1-226
#$ -tc 2
#$ -pe shared 1

module load miniforge/23.11.0
source activate protgps

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_protgps_results
mkdir -p $OUTDIR

TASK=$(printf "%04d" $SGE_TASK_ID)
IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_TSV=$OUTDIR/chunk_${TASK}.tsv

if [ ! -f "$IN_FASTA" ]; then
  echo "No input chunk for task $SGE_TASK_ID, skipping."
  exit 0
fi

START=$(date +%s)
python3 /u/home/d/ddcohn/vep_workflow/localization_condensate/protgps_predict_scale.py $TASK $IN_FASTA $OUT_TSV
RC=$?
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(wc -l < $OUT_TSV 2>/dev/null || echo 0)
echo "Task $SGE_TASK_ID (chunk $TASK): rc=$RC elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$RC" -ne 0 ] || [ "$n_out" -lt "$n_in" ]; then
  exit 1
fi
