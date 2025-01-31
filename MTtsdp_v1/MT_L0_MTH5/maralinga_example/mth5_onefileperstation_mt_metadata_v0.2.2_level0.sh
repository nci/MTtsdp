#!/bin/bash

#PBS -N mth5_full_mt_metadata_v0.2.2_level0
#PBS -q normal
#PBS -P abc
#PBS -l walltime=0:05:00
#PBS -l ncpus=6
#PBS -l mem=24GB
#PBS -l jobfs=1GB
#PBS -l storage=scratch/abc+gdata/up99

module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
mpirun -np $PBS_NCPUS python3.10 ./mth5_onefileperstation_mt_metadata_v0.2.2_level0.py > ./pbs_job_logs/$PBS_JOBID.log

date
