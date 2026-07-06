#!/bin/bash
# Usage: get_last_job.sh [n]
# Prints the job ID(s) of the n most recent SLURM array jobs in logs/.
# n defaults to 1. Output is suitable for use in: sacct -j $(./get_last_job.sh)

n=${1:-1}

jobids=$(ls -t logs/*.out 2>/dev/null \
         | sed 's/.*_\([0-9][0-9]*\)_[0-9][0-9]*\.out$/\1/' \
         | awk '!seen[$0]++' \
         | head -n "$n")

if [ -z "$jobids" ]; then
    echo "No log files found in logs/" >&2
    exit 1
fi

echo "$jobids" | tr '\n' ',' | sed 's/,$/\n/'
