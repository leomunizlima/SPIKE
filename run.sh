#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
LD_LIBRARY_PATH=${HOME}/SPIKE/PARALLEL/pardiso:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
OMP_NUM_THREADS=$9
export OMP_NUM_THREADS
MP_SHARED_MEMORY=yes
export MP_SHARED_MEMORY
mpirun -np $7 ./program $1 1e-8 30 1000 $2 $3 $4 $5 $6 $7 $8 r >trash


