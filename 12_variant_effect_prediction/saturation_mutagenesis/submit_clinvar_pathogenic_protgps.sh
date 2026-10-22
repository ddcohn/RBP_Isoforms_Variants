#!/bin/bash
#$ -N clinvar_pathogenic_protgps
#$ -cwd
#$ -o /u/home/d/ddcohn/clinvar_pathogenic_protgps.log
#$ -j y
#$ -l h_data=12G,h_rt=2:00:00
#$ -pe shared 1

module load miniforge/23.11.0
python3 /u/project/kappel/ddcohn/RBP_pipeline_scripts/saturation_mutagenesis/clinvar_pathogenic_protgps_effects.py
