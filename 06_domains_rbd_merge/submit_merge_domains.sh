#!/bin/bash
#$ -N merge_domains
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_merge_domains.log
#$ -j y
#$ -l h_data=16G,h_rt=1:00:00
#$ -pe shared 1

python3 /u/home/d/ddcohn/_claude_merge_domains_rbd.py
