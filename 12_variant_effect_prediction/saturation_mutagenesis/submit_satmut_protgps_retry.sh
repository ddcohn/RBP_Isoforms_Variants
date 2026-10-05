#!/bin/bash -l
#$ -N protgps_satmut_retry
#$ -cwd
#$ -o /u/home/d/ddcohn/protgps_task_logs/satmut_retry.log.$TASK_ID
#$ -j y
#$ -q kappel_gpu.q
#$ -l h_data=24G,h_rt=6:00:00,highp,require_gpu=1,cuda=1
#$ -t 1-5
#$ -pe shared 1

# Edit CHUNKS below to the task numbers that need redoing. protgps_predict_scale.py
# opens its output file in "w" mode, so a clean resubmit safely overwrites a
# partial/truncated file -- no append or dedup logic needed, unlike the DeepLoc
# retry case. #$ -t above must match CHUNKS' length.
module load miniforge/23.11.0
source activate protgps

CHUNKS=(0077 0113 0147 0174 0183)
TASK=${CHUNKS[$((SGE_TASK_ID - 1))]}

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_protgps_results

IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_TSV=$OUTDIR/chunk_${TASK}.tsv

START=$(date +%s)
python3 /u/home/d/ddcohn/vep_workflow/localization_condensate/protgps_predict_scale.py $TASK $IN_FASTA $OUT_TSV
RC=$?
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(wc -l < $OUT_TSV 2>/dev/null || echo 0)
echo "Retry task $SGE_TASK_ID (chunk $TASK): rc=$RC elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$RC" -ne 0 ] || [ "$n_out" -lt "$n_in" ]; then
  exit 1
fi
