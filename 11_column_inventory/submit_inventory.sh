#!/bin/bash
#$ -N full_inventory
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_inventory.log
#$ -j y
#$ -l h_data=32G,h_rt=2:00:00
#$ -pe shared 1

python3 /u/home/d/ddcohn/_claude_full_inventory.py
