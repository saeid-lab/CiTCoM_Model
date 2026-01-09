#!/bin/bash

# --- Initial setup: create input/output files if missing ---
if [ ! -f all_us.txt ]; then
  echo "Generating all_us.txt..."
  find . -type f -name "us.sub" | sort > all_us.txt
fi

if [ ! -f submitted_us.txt ]; then
  echo "Creating submitted_us.txt..."
  touch submitted_us.txt
fi


# --- Main loop: check queue and submit jobs every 30 minutes ---
#while true; do
while [ -s all_us.txt ]; do
  echo "Checking SLURM queue at $(date)"

  # Count current running or pending jobs
  CURRENT_JOBS=$(squeue -u "$USER" | grep -c " R\| PD")
  echo "Current jobs in queue: $CURRENT_JOBS"

  if [ "$CURRENT_JOBS" -lt 50 ]; then
    echo "Submitting new jobs..."

    # Read all remaining job paths
    mapfile -t JOBS < all_us.txt
    SLICE=("${JOBS[@]:0:250}")
    i=0

    while [ $i -lt ${#SLICE[@]} ]; do
      JOB_PATH="${SLICE[$i]}"
      JOB_DIR=$(dirname "$JOB_PATH")

      if [ -d "$JOB_DIR" ]; then
        SUBMIT_OUTPUT=$(cd "$JOB_DIR" && sbatch us.sub)
        JOB_ID=$(echo "$SUBMIT_OUTPUT" | awk '{print $NF}')

        echo "Job path: $JOB_PATH"
        if [[ -n "$JOB_ID" && "$JOB_ID" =~ ^[0-9]+$ ]]; then
          echo "Submitted job with ID: $JOB_ID"
        else
          echo "Failed to submit job for: $JOB_PATH"
        fi

        # Record the submitted job and remove it from the list
        echo "$JOB_PATH" >> submitted_us.txt
        # grep glitches on the last line, instead using awk!
        # grep -Fxv "$JOB_PATH" all_us.txt > steered.tmp && mv steered.tmp all_us.txt
        awk -v line="$JOB_PATH" '$0 != line' all_us.txt > steered.tmp && mv steered.tmp all_us.txt
        sleep 0.1
      else
        echo "Directory does not exist for path: $JOB_PATH"
      fi

      ((i++))
    done
  else
    echo "Queue full, skipping submission."
  fi

  echo "Sleeping for 15 minutes..."
  sleep 900
done

