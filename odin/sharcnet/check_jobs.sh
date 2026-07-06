#!/bin/bash
# Usage: check_jobs.sh [n]
# Finds the n most recently modified job arrays in logs/ and summarises
# their task states via sacct.  n defaults to 1.

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
    sacct -j "$jobid" --format=State --noheader | sort | uniq -c
done
