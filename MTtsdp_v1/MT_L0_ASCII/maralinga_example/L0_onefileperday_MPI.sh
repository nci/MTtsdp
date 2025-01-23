#!/bin/bash

#PBS -N Maralinga_L0_onefileperday_ASCII
#PBS -q normal
#PBS -P <insert_project_code>
#PBS -l walltime=0:10:00
#PBS -l ncpus=6
#PBS -l mem=24GB
#PBS -l jobfs=10GB
#PBS -l storage=scratch/<insert_project_code>+gdata/up99

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

