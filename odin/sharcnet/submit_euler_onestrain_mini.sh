#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_onestrain_mini
#SBATCH --array=1-24
#SBATCH --time=0-00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/euler_onestrain_mini_%A_%a.out

# Submit from the metapop/odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_onestrain_run_array.R --mini
