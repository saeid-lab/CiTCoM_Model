#!/bin/bash

# summarize_all_sims.sh
# Iterates through simulation folders and prints the latest progress for each.

# Print Header
printf "%-20s %-25s %-8s %-15s %-15s\n" "Directory" "Latest Log" "Part" "Step" "Time"
printf "%-20s %-25s %-8s %-15s %-15s\n" "---------" "----------" "----" "------------" "---------------"

# Iterate over all directories in the current folder
for dir in */; do
    # Remove trailing slash
    dir_name=${dir%/}
    
    # Path to the md folder
    md_path="${dir_name}/md"

    if [ -d "$md_path" ]; then
        # Find the latest log file in that folder using version sort
        LATEST_FILE=$(ls -v "$md_path"/md-*.log 2>/dev/null | tail -n 1)

        if [ -n "$LATEST_FILE" ]; then
            # Extract just the filename for display
            filename=$(basename "$LATEST_FILE")
            
            # Extract Part Number
            part=$(echo "$filename" | grep -oP 'part\K\d+' || echo "N/A")

            # Extract Latest Step and Time
            latest_data=$(tac "$LATEST_FILE" | grep -m 1 -B 1 -E "Step[[:space:]]+Time" | head -n 1)
            step=$(echo "$latest_data" | awk '{print $1}')
            time=$(echo "$latest_data" | awk '{print $2}')

            # Default if no steps recorded yet
            if [[ -z "$step" || "$step" == "Step" ]]; then
                step="--"
                time="--"
            fi

            printf "%-20s %-25s %-8s %-15s %-15s\n" "$dir_name" "$filename" "$part" "$step" "$time"
        else
            # Folder exists but no log files yet
            printf "%-20s %-25s %-8s %-15s %-15s\n" "$dir_name" "No logs found" "--" "--" "--"
        fi
    fi
done
