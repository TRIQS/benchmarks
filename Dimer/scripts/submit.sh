#!/bin/bash
#SBATCH --job-name=benchmark_updates
#SBATCH --nodes=1 --exclusive
#SBATCH --tasks-per-node=128
#SBATCH --constraint=rome
#SBATCH --partition=ccq
#SBATCH --mail-user=nwentzell@flatironinstitute.org
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j

module purge

export MODULEPATH=/mnt/home/wentzell/opt/modules:$MODULEPATH
module load devenv3/clang-py3-mkl

export OMP_NUM_THREADS=1
source $HOME/opt/triqs/share/triqs/triqsvars.sh
export PYTHONPATH=$HOME/opt/pyed:$PYTHONPATH

mpirun cthyb
#mpirun ctint 
#mpirun ctseg
#mpirun w2dyn_cthyb
#mpirun w2dyn_cthyb_delta_interface
