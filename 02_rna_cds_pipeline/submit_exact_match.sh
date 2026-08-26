#!/bin/bash
#$ -N exact_match_rebuild
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_exact_match.log
#$ -j y
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 1

python3 /u/home/d/ddcohn/_claude_exact_match_rebuild.py
