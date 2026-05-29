#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=meta_twostrain
#SBATCH --time=0-02:00:00   # mini2 sacct: ~3 min/task × 14.4 tasks/metajob ≈ 43 min; 2 h gives ~2.8× buffer
#SBATCH --mem=12G            # mini2 MaxRSS ~8.35 GB; 12 G gives ~1.4× buffer
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/meta_twostrain_%j.out

module load meta-farm
module load r/4.5.0
task.run
