#!/bin/bash

nrep=24
for ((i=0; i<nrep; i++)); do
    echo "Entering folder $i"
    cd "$i"
    sbatch eq_rest2.sh
    cd ..
done
