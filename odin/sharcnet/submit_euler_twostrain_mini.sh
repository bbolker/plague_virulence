#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_twostrain_mini
#SBATCH --array=1-96
#SBATCH --time=0-00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/euler_twostrain_mini_%A_%a.out

# Submit from the odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_twostrain_run_array.R --mini
