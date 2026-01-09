## This the prep folder and contains necessary files to prepare Umbrella Sampling MD.
### Required files:
- md_eq7.gro
- steered.mdp
- steered2.mdp
- plumed_steered.dat
- plumed_steered.dat
- job_eq_steered.sh
- job_eq_steered2.sh


### First we need the gro file (md_eq7.gro, md_start.gro). Normally, it comes from previous step after all the prep jobs are done. In here "md_new_nucl2.gro" goes for the next step.
### In order to run job_eq_steered.sh, we need steered.mdp and plumed_steered.dat. In the following you'll see how to prepare the plumed_steered.dat:
Example of 63:
d1: DISTANCE ATOMS=(63A) O2' , (1191M7) C2
d2: DISTANCE ATOMS=(63A) O2' , (1191M7) N2
a1: ANGLE ATOMS=(1191M7) N2, (1191M7) C2, (63A) O2'
a2: ANGLE ATOMS=(1191M7) C2, (63A) C4',   (63A) C2'

Preparing for 64:
awk '$1 ~ /^64[AUCG]?$/' md_eq7.gro | grep "O2'"	=> 2100
awk '$1 ~ /^119/' md_eq7.gro | grep "C2"		=> 3852
awk '$1 ~ /^119/' md_eq7.gro | grep "N2"		=> 3864
awk '$1 ~ /^64[AUCG]?$/' md_eq7.gro | grep "C4'"	=> 2077
awk '$1 ~ /^64[AUCG]?$/' md_eq7.gro | grep "C2'"	=> 2098

d1: DISTANCE ATOMS=2100 , 3852
d2: DISTANCE ATOMS=2100 , 3864
a1: ANGLE ATOMS=3864, 3852, 2100
a2: ANGLE ATOMS=3852, 2077, 2098

### Now, making the index.ndx file:
module load  .archive/CCE-GPU-3.0.0 gromacs/2023_amd-mpi-omp-plumed-python3
gmx_mpi make_ndx -f md_eq7.gro # just q
(I copied it from previous step)

### After that: 
sbatch job_eq_steered.sh

### After that:
sbatch job_eq_steered2.sh


## That's all we need for this folder.

