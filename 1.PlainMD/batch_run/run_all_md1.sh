#!/bin/bash

# Load modules once
module load .archive/CCE-GPU-3.0.0 gromacs/2023_amd-mpi-omp-plumed-python3

for dir in */; do
    md_dir="${dir}md"

    # Skip if md directory does not exist
    [ -d "$md_dir" ] || continue

    cd "$md_dir" || continue

    echo 'q' | gmx_mpi make_ndx -f min.gro

    gmx_mpi grompp \
        -f ../mdp_single/md-2.mdp \
        -c md_eq6.gro \
        -p ../top/topol.top \
        -o md-2 \
        -n index \
        -maxwarn 1

    sbatch job_restart.sh

    cd - > /dev/null
done
