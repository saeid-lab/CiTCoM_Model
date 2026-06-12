#!/bin/bash

for dir in */; do
    cd "$dir/md" && sbatch auto_md_eq.sh && cd -
done
