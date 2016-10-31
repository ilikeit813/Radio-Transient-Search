#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
#PBS -W group_list=blueridge
#PBS -q normal_q
#PBS -A LWAradio
#PBS -M ilikeit@vt.edu
#PBS -m bea

module reset
module swap mvapich2 openmpi
module load mkl python

cd /work/blueridge/ilikeit/057139_000656029/
cp /home/ilikeit/ft.py .
cp /home/ilikeit/errors.py .
cp /home/ilikeit/drx.py .
cp /home/ilikeit/dp.py .
#rm *.txt
mpirun -np $PBS_NP python ft.py 057139_000656029 1

echo "------------------------------------------"
echo "done! with cores =" $PBS_NP
echo "------------------------------------------"

exit;
