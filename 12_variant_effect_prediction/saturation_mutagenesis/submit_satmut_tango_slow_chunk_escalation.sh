#!/bin/bash -l
#$ -N tango_satmut_slow
#$ -cwd
#$ -o /u/home/d/ddcohn/tango_task_logs/satmut_slow.log.$TASK_ID
#$ -j y
#$ -q kappel_cpu.q
#$ -l h_data=2G,h_rt=168:00:00,highp
#$ -t 1-3
#$ -pe shared 1

# For chunks that are consistently too slow to finish within the general
# pool's 24h cap even after multiple resubmissions (identified by comparing
# elapsed-time-since-job-start against row-count progress across two checks:
# a rate well under the ~1.4 rows/sec baseline means the chunk is dominated
# by unusually long sequences -- TANGO's cost scales worse than linearly
# with length, up to ~21x slower at the longest RBP in this project's list
# vs. a typical-length one). Escalates to kappel_cpu.q, which allows up to
# 336h (14 days) -- set h_rt here to whatever headroom the observed rate
# implies, this used 168h (7 days) as a safe margin over an estimated ~38-59h
# remaining need. tango_predict_chunk.py is resumable at the sequence level,
# so this continues exactly where a prior timed-out attempt left off --
# no progress lost, same mechanism as any other retry.
#
# Edit CHUNKS below to the task numbers that need escalating. #$ -t above
# must match CHUNKS' length.
#
# Don't try to extend an already-running task's h_rt live via qalter --
# tested, doesn't work. `qalter -l h_data=2G,h_rt=HH:MM:SS,highp <job_id>`
# reports success and qstat -j immediately shows the new h_rt, but the
# execution daemon's actual wall-time enforcement on the already-dispatched
# task keeps using the ORIGINAL limit -- confirmed via qacct showing a kill
# at exactly the old h_rt with "execd enforced h_rt limit", despite the
# scheduler metadata claiming the job had 7 days left. Just resubmit this
# script fresh with the real h_rt set from the start; resumability means
# nothing is lost beyond whatever time the failed live-extension wasted.
module load miniforge/23.11.0

CHUNKS=(0184 0188 0219)
TASK=${CHUNKS[$((SGE_TASK_ID - 1))]}

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_tango_results
TANGO_BIN=/u/home/d/ddcohn/tango_x86_64_release
DRIVER=/u/home/d/ddcohn/_claude_tango_predict_chunk.py

IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_TSV=$OUTDIR/chunk_${TASK}.tsv

START=$(date +%s)
python3 $DRIVER $IN_FASTA $OUT_TSV $TANGO_BIN
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(($(wc -l < $OUT_TSV 2>/dev/null || echo 0) - 1))
echo "Slow-chunk task $SGE_TASK_ID (chunk $TASK): elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$n_out" -lt "$n_in" ]; then exit 1; fi
