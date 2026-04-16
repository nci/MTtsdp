#!/bin/bash

#PBS -N North_Flinders_L0_entire_ASCII
#PBS -q normal
#PBS -P <aa0>
#PBS -l walltime=0:10:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/up99+gdata/<xy12>
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
mpirun -np $PBS_NCPUS python3.12 ./02_merge_singlestation.py

date
