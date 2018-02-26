#!/bin/bash
# CCY - Added OPT-INSTALL_DIR as a template variable that is replaced by my git-export.sh script with
# the path to which the package is installed.  If you don't have my git-export.sh script, then you can
# just manually replace all instances of OPT-INSTALL_DIR as necessary.
#

USAGE='waterfall.sh

   waterfall.sh [--help | -h] <data_path>

   Script calling the waterfall.py python script in the installed package of radiotrans to operate on a
   radio time series specified by <data_path>.  The results of waterfall.py are output to the current
   working directory.

   OPTIONS:
   --help | -h : Display this help message.
  ' 

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

#module reset
#module swap openmpi
#module swap mvapich2 openmpi
#module load mkl python openmpi
NUM_PROCS=0



# Path to the data file.
DATA_PATH=

# Get user input from command-line.
while [ -z "${DATA_PATH}" -a -n "$1" ]
do
   case  "$1"
   in
      --help | - h) # Display usage and then exit
         echo "${USAGE}"
         exit 0
      -*) # Any other option flag is ignored
         echo "UNKNOWN OPTION: $1"
         echo "${USAGE}"
         exit 1
      *) # Use the first non-option flag as the data path.
         DATA_PATH="$1"
         shift
         ;;
   esac
done

if [ -n "${DATA_PATH}" ]; then
   echo "running waterfall.py with DATA_PATH=${DATA_PATH}"
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/waterfall.py "${DATA_PATH}"
   echo "done"
else
   echo "No data path specified"
   echo "${USAGE}"
   exit 1
fi

exit 0;
