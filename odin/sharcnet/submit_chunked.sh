#!/bin/bash
# Submit a large job array in chunks to stay under the account's
# AssocMaxSubmitJobLimit (pending+running jobs per user, counted across
# all array tasks). Waits for queue headroom before each chunk.
#
# Usage: submit_chunked.sh <submit_script> <first> <last> [chunk_size] [max_in_queue] [poll_seconds]
#
# Example (8100-task array, 1000-job account cap):
#   bash submit_chunked.sh submit_euler_twostrain_singlepatchintro.sh 1 8100
#
# Run under tmux/screen or nohup — this can take hours to drain all chunks.

submit_script=$1
first=$2
last=$3
chunk_size=${4:-900}
max_in_queue=${5:-950}
poll_seconds=${6:-300}

if [ -z "$submit_script" ] || [ -z "$first" ] || [ -z "$last" ]; then
    echo "Usage: $0 <submit_script> <first> <last> [chunk_size] [max_in_queue] [poll_seconds]" >&2
    exit 1
fi

wait_for_room () {
    needed=$1
    while true; do
        n_queued=$(squeue -u "$USER" -h -r -t pending,running | wc -l)
        if [ "$((n_queued + needed))" -le "$max_in_queue" ]; then
            return
        fi
        echo "$(date '+%F %T'): ${n_queued} jobs already queued/running; waiting for room for ${needed} more..."
        sleep "$poll_seconds"
    done
}

start=$first
while [ "$start" -le "$last" ]; do
    end=$((start + chunk_size - 1))
    if [ "$end" -gt "$last" ]; then
        end=$last
    fi
    chunk_n=$((end - start + 1))

    wait_for_room "$chunk_n"

    echo "$(date '+%F %T'): submitting --array=${start}-${end} (${chunk_n} tasks)"
    if sbatch --array="${start}-${end}" "$submit_script"; then
        start=$((end + 1))
    else
        echo "$(date '+%F %T'): sbatch failed for --array=${start}-${end}; will retry after a pause"
        sleep "$poll_seconds"
    fi
done

echo "$(date '+%F %T'): all chunks submitted (${first}-${last})."
