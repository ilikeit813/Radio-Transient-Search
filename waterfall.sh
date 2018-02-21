#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=2:ppn=6
#PBS -W group_list=hokieone
#PBS -q normal_q
#PBS -A hokieone

# Add any the intel compiler and MPT MPI modules
#module reset
#module load mkl mpiblast python
#module swap intel gcc
#module load mkl python
#module load mkl mpt python
#module load intel mpt
#module add mpiblast 
#module load openmpi

#module load mkl mpiblast python

module reset
#module swap openmpi
#module swap mvapich2 openmpi
module load mkl python openmpi
NUM_PROCS=1

DATA_DIR="/data/network/recent_data/jtsai"
DATA_FILENAME="057974_001488582"
DATA_PATH="${DATA_DIR}/${DATA_FILENAME}"

mpirun -np ${NUM_PROCS} python waterfall.py "${DATA_PATH}"

echo "done"
exit;
