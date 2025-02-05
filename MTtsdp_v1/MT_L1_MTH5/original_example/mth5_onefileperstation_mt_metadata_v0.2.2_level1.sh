#!/bin/bash

#PBS -N mth5_full_mt_metadata_v0.2.2_level1
#PBS -q normal
#PBS -P <insert_compute_project_code>
#PBS -l walltime=0:05:00
#PBS -l ncpus=96
#PBS -l mem=384GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/fp0+gdata/my80+gdata/lm70+gdata/up99

### current stable version of NCI-geophys module
#module use /g/data/up99/modulefiles
#module load NCI-geophys/22.11

### latest version of NCI-geophys with mt_metadata v0.2.2. Not yet released to up99.
module use /g/data/fp0/apps/Modules/modulefiles
module load NCI-geophys/23.02

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
mpirun -np $PBS_NCPUS python3 ./mth5_onefileperstation_mt_metadata_v0.2.2_level1.py > ./pbs_job_logs/$PBS_JOBID.log

date

