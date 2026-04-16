#!/bin/bash

#PBS -N NE_SA_L0_ASCII
#PBS -q normal
#PBS -P <abc>
#PBS -l walltime=0:15:00
#PBS -l ncpus=12
#PBS -l mem=48GB
#PBS -l jobfs=1GB
#PBS -l storage=gdata/<xy12>+gdata/up99
#PBS -j oe
#PBS -o /g/data/<xy12>/supplementing_files/pbs_job_logs

python_script="/g/data/<xy12>/<abc123>/used_scripts/MT_L0_ASCII/L0_onefileperday/L0_onefileperday_MPI.py"
log_dir="/g/data/<xy12>/<abc123>/supplementing_files/pbs_job_logs"

# mkdir -p ${log_dir}

### Load NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/25.03

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
mpirun -np $PBS_NCPUS python3.12 ${python_script} > ${log_dir}/$PBS_JOBID.log

date

