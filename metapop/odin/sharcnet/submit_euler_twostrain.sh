#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_twostrain
#SBATCH --array=1-56000
#SBATCH --time=0-00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/euler_twostrain_%A_%a.out

# Submit from the metapop/odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs
# NOTE: 56000 tasks may exceed the cluster's MaxArraySize. If so, submit
# in slices: --array=1-1000, --array=1001-2000, etc.

module load r/4.5.0

Rscript euler_twostrain_run_array.R
