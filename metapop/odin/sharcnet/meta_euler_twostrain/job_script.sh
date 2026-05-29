#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=meta_twostrain
#SBATCH --time=0-08:00:00   # est. ~14 min/task (scaled from mini: 6.5 s × 4 CPU / 50 sims × 2 dt × 4 patches × 200 sims) → ~3.4 h/metajob worst case; 8 h adds ~2.4× buffer
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/meta_twostrain_%j.out

module load meta-farm
module load r/4.5.0
task.run
