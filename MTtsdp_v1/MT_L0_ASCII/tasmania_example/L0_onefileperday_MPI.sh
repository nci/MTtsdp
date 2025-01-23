#!/bin/bash

#PBS -N Tas_L0_onefileperday_ASCII
#PBS -q normal
#PBS -P abc
#PBS -l walltime=0:10:00
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/abc+gdata/my80+gdata/up99+scratch/abc

### current stable version of NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
mpirun -np $PBS_NCPUS python3.10 ./L0_onefileperday_MPI.py > ./pbs_job_logs/$PBS_JOBID.log

date

