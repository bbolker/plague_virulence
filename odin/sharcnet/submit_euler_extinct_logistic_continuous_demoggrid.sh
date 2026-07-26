#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=extinct_logcon_dg
#SBATCH --array=1-2600
#SBATCH --time=0-00:05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/euler_%A_%a.out

# Submit from the odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs
# 2600 tasks exceeds the account's MaxSubmitJobsPerUser (1000 on Nibi) —
# submit via submit_chunked.sh instead of sbatch directly:
#   bash submit_chunked.sh submit_euler_extinct_logistic_continuous_demoggrid.sh 1 2600

module load r/4.5.0

Rscript euler_onepatch_onestrain_extinct_run_array.R --demog_grid
