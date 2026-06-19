#!/bin/bash

#PBS -N mth5_full_mt_metadata_v0.2.2_level1
#PBS -q normalsr
#PBS -P abc
#PBS -l walltime=0:05:00
#PBS -l ncpus=30
#PBS -l mem=120GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/ab12+gdata/up99

module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo

SURVEY="SE_SA_survey"
STATIONS_SINGLE_RUN="['SA01','SA02','SA03','SA03_2','SA04','SA05','SA06','SA07','SA08','SA09','SA09_2','SA10','SA13','SA14','SA15','SA19_1','SA19_2','SA20','SA22','SA25']"
STATIONS_MULTIPLE_RUNS="['SA11','SA12','SA16','SA17','SA18','SA19_3','SA21','SA23']"


mpirun -np $PBS_NCPUS python3 ./mth5_onefileperstation_mt_metadata_v0.2.2_level1.py "$SURVEY" "$STATIONS_SINGLE_RUN" "$STATIONS_MULTIPLE_RUNS" > ./pbs_job_logs/$PBS_JOBID.log

date

