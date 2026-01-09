#!/bin/bash

# --- Initial setup: create input/output files if missing ---
if [ ! -f all_steered.txt ]; then
  echo "Generating all_steered.txt..."
  find . -type f -name "steered.sub" | sort > all_steered.txt
fi

if [ ! -f submitted_steered.txt ]; then
  echo "Creating submitted_steered.txt..."
  touch submitted_steered.txt
fi


# --- Main loop: check queue and submit jobs every 30 minutes ---
while true; do
  echo "Checking SLURM queue at $(date)"

  # Count current running or pending jobs
  CURRENT_JOBS=$(squeue -u "$USER" | grep -c " R\| PD")
  echo "Current jobs in queue: $CURRENT_JOBS"

  if [ "$CURRENT_JOBS" -lt 50 ]; then
    echo "Submitting new jobs..."

    # Read all remaining job paths
    mapfile -t JOBS < all_steered.txt
    SLICE=("${JOBS[@]:0:250}")
    i=0

    while [ $i -lt ${#SLICE[@]} ]; do
      JOB_PATH="${SLICE[$i]}"
      JOB_DIR=$(dirname "$JOB_PATH")

      if [ -d "$JOB_DIR" ]; then
        SUBMIT_OUTPUT=$(cd "$JOB_DIR" && sbatch steered.sub)
        JOB_ID=$(echo "$SUBMIT_OUTPUT" | awk '{print $NF}')

        echo "Job path: $JOB_PATH"
        if [[ -n "$JOB_ID" && "$JOB_ID" =~ ^[0-9]+$ ]]; then
          echo "Submitted job with ID: $JOB_ID"
        else
          echo "Failed to submit job for: $JOB_PATH"
        fi

        # Record the submitted job and remove it from the list
        echo "$JOB_PATH" >> submitted_steered.txt
        # grep glitches on the last line, instead using awk!
        # grep -Fxv "$JOB_PATH" all_steered.txt > steered.tmp && mv steered.tmp all_steered.txt
        awk -v line="$JOB_PATH" '$0 != line' all_steered.txt > steered.tmp && mv steered.tmp all_steered.txt
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

