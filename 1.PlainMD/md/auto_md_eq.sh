#!/bin/bash
#SBATCH --account=ibc7585
#SBATCH --job-name=md_eq_*
#SBATCH --constraint=MI250
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --output=cines_%x_%J.out
#SBATCH --error=cines_%x_%J.err
#SBATCH --time=24:00:00


module purge

nodes="${SLURM_JOB_NUM_NODES}"
ppn="${SLURM_NTASKS_PER_NODE}"
threads="${BM_THREADS:-$SLURM_CPUS_PER_TASK}"
smt=1
#ranks=$(( nodes * ppn ))
ranks="${SLURM_NTASKS}"

system=adastra
cpu=7A53
arch=trento
gpu=mi250
build=cce
echo "System: $system, Arch: $arch, CPU: $cpu, GPU: $gpu"
echo "Ranks: $ranks, PPN: $ppn, Nodes: $nodes, Threads: $threads, SMT: $smt, GPN: $gpn"
echo "Build: $build"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on: $SLURM_JOB_NODELIST"

# Stack limits
ulimit -c unlimited
ulimit -s unlimited
export OMP_NUM_THREADS="$threads"
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Gromacs envars
module load  .archive/CCE-GPU-3.0.0
module load gromacs/2023_amd-mpi-omp-plumed-python3


export GMX_ENABLE_DIRECT_GPU_COMM=yes
export GMX_DISABLE_ALTERNATING_GPU_WAIT=yes
export GMX_ENABLE_STAGED_GPU_TO_CPU_PMEPP_COMM=yes
#
##
export GMX_GPU_DD_COMMS=1
export GMX_GPU_PME_PP_COMMS=1
export GMX_NO_QUOTES=1

export OMP_NUM_THREADS=$threads
export SLURM_CPU_BIND_TYPE=rank
export CRAY_ACC_DEBUG=3
env OMP_NUM_THREADS=$threads
#############################
export GMX_FORCE_GPU_AWARE_MPI=yes
export MPICH_OFI_NIC_POLICY=GPU
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export GMX_DISABLE_GPU_TIMING=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
#############################


# Define variables
MDP_PATH="./../mdp_single"
TOP_PATH="./../top"
INDEX_FILE="$TOP_PATH/index.ndx"
TOPOLOGY="$TOP_PATH/topol.top"

# Define the sequence of equilibration steps
declare -a mdp_files=("relax_wat.mdp" "heat.mdp" "md_eq1.mdp" "md_eq2.mdp" "md_eq3.mdp" "md_eq4.mdp" "md_eq5.mdp" "md_eq6.mdp")
declare -a gro_files=("min.gro" "relax_wat.gro" "heat.gro" "md_eq1.gro" "md_eq2.gro" "md_eq3.gro" "md_eq4.gro" "md_eq5.gro")
declare -a output_files=("relax_wat" "heat" "md_eq1" "md_eq2" "md_eq3" "md_eq4" "md_eq5" "md_eq6")

# Loop through equilibration steps
for i in {0..7}; do
    filein="$MDP_PATH/${mdp_files[$i]}"
    filegro="${gro_files[$i]}"
    file="${output_files[$i]}"

    echo "Running equilibration step: $file"
    gmx_mpi grompp -f $filein -c $filegro -p $TOPOLOGY -o $file -n $INDEX_FILE -maxwarn 1 -r $filegro
    
    srun -K1 --mpi=cray_shasta -N $SLURM_NNODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
        --accel-bind verbose --cpu-bind cores -m block:block \
        gmx_mpi mdrun -s $file.tpr -v -deffnm $file -ntomp $SLURM_CPUS_PER_TASK \
        -pin on -nb gpu -bonded gpu -update cpu -maxh 23

done

echo "Equilibration completed successfully."

