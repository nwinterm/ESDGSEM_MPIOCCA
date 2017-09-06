#!/bin/bash -l
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=32
#SBATCH --mem=16gb
#SBATCH --time=80:00:00
#SBATCH --account=AG-Gass

# number of nodes in $SLURM_NNODES (default: 1)
# number of tasks in $SLURM_NTASKS (default: 1)
# number of tasks per node in $SLURM_NTASKS_PER_NODE (default: 1)
# number of threads per task in $SLURM_CPUS_PER_TASK (default: 1)

module load openmpi
export OCCA_DIR=/home/nwinter2/occa-master
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib
export PATH=$PATH:$OCCA_DIR/bin
srun -n $SLURM_NTASKS main Serial >& output.txt

