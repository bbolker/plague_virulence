#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=extinct_logrf
#SBATCH --array=1-520
#SBATCH --time=0-00:05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/euler_%A_%a.out

# Submit from the odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_onepatch_onestrain_extinct_run_array.R --reedfrost
