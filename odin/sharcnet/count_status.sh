#!/bin/bash
# Usage: count_status.sh [n]
# Like check_jobs.sh but counts one record per array task (filters out .batch sub-steps).

n=${1:-1}

jobids=$(ls -t logs/*.out 2>/dev/null \
         | sed 's/.*_\([0-9][0-9]*\)_[0-9][0-9]*\.out$/\1/' \
         | awk '!seen[$0]++' \
         | head -n "$n")

if [ -z "$jobids" ]; then
    echo "No log files found in logs/" >&2
    exit 1
fi

for jobid in $jobids; do
    echo "=== Job $jobid ==="
    sacct -j "$jobid" --format=JobID%30,State --noheader \
        | awk '$1 ~ /_[0-9]+$/ {print $2}' \
        | sort | uniq -c
done
