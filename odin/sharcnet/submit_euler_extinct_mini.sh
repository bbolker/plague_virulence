#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_extinct_mini
#SBATCH --array=1-32
#SBATCH --time=0-00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/euler_mini_%A_%a.out

# Submit from the odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_onepatch_onestrain_extinct_run_array_mini.R
