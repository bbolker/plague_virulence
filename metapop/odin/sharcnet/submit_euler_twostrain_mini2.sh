#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=twostrain_mini2
#SBATCH --array=1-900
#SBATCH --time=0-00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/euler_twostrain_mini2_%A_%a.out

# Submit from the metapop/odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_twostrain_run_array.R --mini2
