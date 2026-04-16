#!/bin/bash

#PBS -N Maralinga_LEMI_L0_ASCII_perrun
#PBS -q normal
#PBS -P <aa0>
#PBS -l walltime=0:10:00
#PBS -l ncpus=2
#PBS -l mem=8GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/<xy12>+gdata/up99
#PBS -j oe
#PBS -o /g/data/<xy12>/<abc123>/supplementing_files/pbs_job_logs

### Load NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/25.03

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
mpirun -np $PBS_NCPUS python3.12 ./02_merge.py

date
