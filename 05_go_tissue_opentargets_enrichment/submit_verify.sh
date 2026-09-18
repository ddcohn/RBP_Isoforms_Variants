#!/bin/bash
#$ -N verify_merge
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_verify.log
#$ -j y
#$ -l h_data=16G,h_rt=1:00:00
#$ -pe shared 1

python3 /u/project/kappel/ddcohn/RBP/scripts/verify_merge.py
