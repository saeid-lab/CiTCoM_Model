#!/bin/bash
#SBATCH --account=ibc7585
#SBATCH --job-name=Production_*
#SBATCH --constraint=MI250
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --output=cines_%x_%J.out
#SBATCH --error=cines_%x_%J.err
#SBATCH --time=24:00:00


ROOTNAME=md-2              # Rootname for gromacs files
MDRUN_OPT=""         # Additional options to use in mdrun
BATCH_FNAME=job_restart.sh   # The file name of the submition file (this file)
CHECK_DURATION=120           # Time to count for the post-MD checks (in seconds)

# You should not have to change anything after this point
# (except if you want to change gromacs version, then you should change how
# gromacs is loaded)

### Load gromacs

#########################################################
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


MDRUN="srun -K1 --mpi=cray_shasta  -N $nodes -n $ranks -c $threads  --accel-bind verbose  --cpu-bind cores -m block:block gmx_mpi mdrun  "


### Read the walltime and convert it in maxh value
# We need to keep some time after the MD run for the checks,
# so we substract CHECK_DURATION from the walltime to get maxh.
echo "jobid is $SLURM_JOB_ID "
echo $(scontrol show jobid $SLURM_JOB_ID)
walltime=$(scontrol show jobid $SLURM_JOB_ID | \
           tr ' ' '\n' | \
           grep "TimeLimit" | \
           cut -f 2 -d = | \
           awk -F '-' -v check_duration=$CHECK_DURATION \
               'BEGIN{seconds=0} \
                { \
                    if (NF==1) {days=0; hms=$1} \
                    else {days=$1; hms=$2}; \
                    seconds=seconds + (days * 24 * 3600); \
                    time_split_length=split(hms, time_split, ":"); \
                    for (i=time_split_length; i>0; i--) { \
                        seconds=seconds + \
                                (time_split[i] \
                                     * 60**(time_split_length - i)) \
                    } \
                } \
                END{print (seconds - check_duration)/3600}')
echo "walltime is $walltime"

### Move in the working directory
# Make sure any symbolic links are resolved to absolute path
export WORKDIR=$(readlink -f $SLURM_SUBMIT_DIR)
echo "Work directory is $WORKDIR"
  
# Change to the direcotry that the job was submitted from
cd $WORKDIR

### Identify the run
# The information on the run number is stored in the last_cycle file.
# This file needs to be created if it does not already exist.
if [[ ! -e last_cycle ]]
then
    echo 0 > last_cycle
fi

# At this point, the last_cycle file exists. We read the index of the last run
# and we increment the number
prev_cycle=$(tail -n1 last_cycle)
cycle=$(( $prev_cycle + 1 ))

echo "The current cycle is $cycle"

# Update the last_cycle file for the next round
echo $cycle >> last_cycle

### Identify the previous run
# By default, the previous run is the one we read in the last_cycle file. If a
# run crashed before it writes a checkpoint, then we cannot start from it. So
# we need to find the last checkpoint available. If there is no checkpoint,
# then prev_cycle is 0 and we need to start from the beginning.
while [[ (! -e "$ROOTNAME.${prev_cycle}.cpt") && (${prev_cycle} -gt 0) ]]
do
    let prev_cycle--
done
checkpoint=$ROOTNAME.$prev_cycle.cpt

# Sometime the simulation started with an other script and the checkpoint is
# not numbered. We still want to continue from this checkpoint.
if [[ ($prev_cycle -eq 0) && (-e $ROOTNAME.cpt) ]]
then 
    checkpoint=$ROOTNAME.cpt
fi

echo "We will use this checkpoint: $checkpoint"

### Is the simulation finished already?
# If the simulation already reached the number of steps requested in the TPR,
# then it is useless to start a new run.
# We first read the TRP file to know how many steps were requested.
tpr_nsteps=$(srun gmx_mpi dump -s ${ROOTNAME}.tpr 2> /dev/null | grep nsteps | cut -f 2 -d = | sed 's/[^0-9]//')
#tpr_nsteps=$(gmx_mpi dump -s ${ROOTNAME}.tpr 2> /dev/null | \            grep nsteps | cut -f 2 -d = | sed 's/[^0-9]//')

# We need to find what is the last step that has been simulated. We read it
# from the log file of the previous cycle.
last_step=$(grep "Writing checkpoint" $ROOTNAME.$prev_cycle.part*.log | tail -n1 | cut -f 4 -d ' ')
echo "Requested number of steps: $tpr_nsteps"
echo "Last step of the previous run: $last_step"

### Run the simulation if appropriate
# We do the run if we are not done or if there is no previous log file. This
# last case can happen if (1) it is the first run or (2) we start from a run
# done with an other script. If we start from a run that used an other script,
# then we try anyway because it is painful to detect and it will not cost much
# anyway.
if [[ ($(ls $ROOTNAME.$prev_cycle.part*.log | wc -l) -eq 0) || \
      ($last_step -lt $tpr_nsteps) ]]
then
    # Launch the parallel job
    $MDRUN  -ntomp $threads -pin on -nb gpu -bonded gpu -update cpu  -nice 0 -s $ROOTNAME -deffnm $ROOTNAME.$cycle -v  -stepout 1000 -maxh 23 -cpi $checkpoint -noappend  $MDRUN_OPT >& $ROOTNAME.$cycle.runout
    
    # Check if the trajectory and the energy file are not corrupted. We will
    # not relaunch a run if a corruption occured.
    srun -n 1 gmx_mpi check -f $ROOTNAME.$cycle.part*.xtc
    #gmx_mpi check -f $ROOTNAME.$cycle.part*.xtc
    check_xtc=$?
    echo "XTC check output code is $check_xtc"
    srun -n 1 gmx_mpi check -e $ROOTNAME.$cycle.part*.edr
    #gmx_mpi check -e $ROOTNAME.$cycle.part*.edr
    check_edr=$?
    echo "EDR check output code is $check_edr"
    integrity=$(( $check_xtc + $check_edr ))
    if [[ $integrity -ne 0 ]]
    then
        echo "Error with XTC or EDR. Check for corruption."
    else
        echo "Integrity OK"
    fi
 
    # If we are not done, then we need to requeue a job
    last_step=$(grep "Writing checkpoint" $ROOTNAME.$cycle.part*.log | \
            tail -n1 | cut -f 4 -d ' ')
    echo "Last simulated step is $last_step"
    if [[ $last_step -lt $tpr_nsteps && $integrity -eq 0 ]]
    then
        sbatch $BATCH_FNAME
    fi
fi

echo "Run $cycle is done"
