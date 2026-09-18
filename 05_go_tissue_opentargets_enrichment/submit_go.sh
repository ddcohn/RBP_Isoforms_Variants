#!/bin/bash
#$ -N go_pull
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_go.log
#$ -j y
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 1

python3 /u/project/kappel/ddcohn/RBP/scripts/pull_go.py
