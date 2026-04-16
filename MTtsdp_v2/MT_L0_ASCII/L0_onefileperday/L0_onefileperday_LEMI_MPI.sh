#!/bin/bash

#PBS -N LEMI_onefileperday_ASCII
#PBS -q normal
#PBS -P <abc>
#PBS -l walltime=0:05:00
#PBS -l ncpus=14
#PBS -l mem=56GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/<xy12>+gdata/up99

### NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/25.03

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
mpirun -np $PBS_NCPUS python3.12 ./L0_onefileperday_LEMI_MPI.py > ./logs/$PBS_JOBID.log

date
