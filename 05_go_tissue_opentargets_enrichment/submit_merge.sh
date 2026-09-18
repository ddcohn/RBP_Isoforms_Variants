#!/bin/bash
#$ -N merge_all
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_merge.log
#$ -j y
#$ -l h_data=32G,h_rt=2:00:00
#$ -pe shared 1

python3 /u/project/kappel/ddcohn/RBP/scripts/merge_all_pending.py
