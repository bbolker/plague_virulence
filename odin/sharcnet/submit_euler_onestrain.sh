#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_onestrain
#SBATCH --array=1-1400
#SBATCH --time=0-00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/euler_onestrain_%A_%a.out

# Submit from the metapop/odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_onestrain_run_array.R
