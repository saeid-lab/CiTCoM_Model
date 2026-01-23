#!/usr/bin/env bash

OUT_FOLDER="MLv3-1ns"

source /usr/local/anaconda3/etc/profile.d/conda.sh
conda activate mdtraj

for dir in 1AUD 1EHZ 1GID 1NBS 1VTQ 1Y26 2GDI 2K95 2L1V 2L8H 3DIG 3IVK 3MXH 3P49 3PDR
do
    target="$dir/Results/$OUT_FOLDER"
    if [ ! -d "$target" ]; then
        echo "Creating $target"
        mkdir -p "$target"
    fi
    
    if [ -d "$target" ]; then
        echo "Processing $target"
        (
            cd "$target" || exit 1

            # Run trjconv only if pdbs.pdb does NOT exist
            if [ ! -f pdbs.pdb ]; then
                echo "pdbs.pdb not found, generating..."
                echo 'RNA' | gmx trjconv -s ../md-2.tpr -f "../md-${dir}_whole_fit_Rep123.xtc" -n ../index.ndx -o pdbs.pdb -skip 500
            else
                echo "pdbs.pdb already exists, skipping gmx trjconv"
            fi

            # Run x3dna-dssr only if pdbs.json does NOT exist
            if [ ! -f pdbs.json ]; then
                echo "pdbs.json not found, generating..."
                x3dna-dssr --input=pdbs.pdb --output=pdbs.json --json --prefix=traj --md
            else
                echo "pdbs.json already exists, skipping x3dna-dssr"
            fi

            python /media/saeid/Disk_2/RNAs/SHAPER/V0_FINAL_with_atomic_angles_10.py
            python /media/saeid/Disk_2/RNAs/SHAPER/correlation.py
        )
        echo "Done with $target"
    else
        echo -e "\e[0;31mSkipping $target - does not exist\e[0m"
    fi
done


# Confirming RNAs, one by one, to ensure consistence length between N-nuc and N-Reac AND valid Reac values
# 1AUD 1EHZ 1GID 1Y26 2L8H 3DIG 3IVK
# 1NBS - missing 5 nucs
# 1VTQ - missing 1 nuc
# 2GDI - missing 2 nucs
# 2K95 - missing 1 nucs
# 2L1V - missing 2 nucs
# 3MXH - missing 1 nuc
# 3P49 - missing 10 nucs
# 3PDR - missing 7 nucs
