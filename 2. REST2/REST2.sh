#!/bin/bash
#SBATCH --account=ibc7585
#SBATCH --job-name="REST2_1F84"
#SBATCH --constraint=MI250
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=8         # 1 par GPU
#SBATCH --cpus-per-task=8           # nombre de coeurs ?|  réserver
#SBATCH --gpus-per-node=8
#SBATCH --exclusive         # mandatory with mps
#SBATCH --hint=nomultithread
##SBATCH -C mps              # to allow GPU sharing between processes
#SBATCH --output=cines_%x_%J.out
#SBATCH --error=cines_%x_%J.err
#SBATCH --time=00:20:00

# nettoyage des modules charges en interactif et hérités par défaut
module purge

nodes="${SLURM_JOB_NUM_NODES}"
ppn="${SLURM_NTASKS_PER_NODE}"
threads="${BM_THREADS:-$SLURM_CPUS_PER_TASK}"
smt=1
ranks="${SLURM_NTASKS}"
echo "Ranks: $ranks, PPN: $ppn, Nodes: $nodes, Threads: $threads, SMT: $smt, GPN: $gpn"
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
##
export GMX_GPU_DD_COMMS=1
export GMX_GPU_PME_PP_COMMS=1
export GMX_NO_QUOTES=1
export SLURM_CPU_BIND_TYPE=rank
export CRAY_ACC_DEBUG=3

srun -K1 --mpi=cray_shasta -N $nodes -n $ranks -c $threads --accel-bind verbose --gpu-bind=verbose,closest -m block:block gmx_mpi mdrun -multidir {0..23} -deffnm run -plumed plumed.dat -hrex -replex 1000 -ntomp $threads -pin on -nb gpu -pme gpu -bonded gpu -update cpu -cpi run.cpt -maxh 23
