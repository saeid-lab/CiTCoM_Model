#!/bin/bash

# === Configuration ===
start_line=1938
end_line=2907
input_file="processed.top"
output_file="modified.top"

# === Processing ===
awk -v start="$start_line" -v end="$end_line" '
function pad(val, len) {
    # Pad or trim to exactly "len" characters
    fmt = "%-" len "s"
    return sprintf(fmt, val)
}
{
    if (NR >= start && NR <= end && $0 !~ /^;/ && NF >= 2) {
        # Capture all fields
        orig_line = $0
        field1 = $1
        atom = $2

        # Add _ only if not already there
        if (atom !~ /_$/) {
            atom = atom "_"
        }

        # Rebuild line with exact field widths (align spacing with padding)
        # Adjust widths based on visual inspection; tweak if needed
        printf "%6s%10s", field1, atom
        for (i = 3; i <= NF; i++) {
            printf "%8s", $i
        }
        printf "\n"
    } else {
        print
    }
}' "$input_file" > "$output_file"

rm $input_file
mv modified.top processed.top

