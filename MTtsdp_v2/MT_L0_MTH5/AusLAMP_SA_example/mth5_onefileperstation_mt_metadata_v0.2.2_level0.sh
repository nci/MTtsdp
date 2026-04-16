#!/bin/bash

#PBS -N mth5_full_mt_metadata_v0.2.2_level0_LEMI
#PBS -q normalsr
#PBS -P <a00>
#PBS -l walltime=0:10:00
#PBS -l ncpus=4
#PBS -l mem=40GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/<xy12>+gdata/up99

module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo

SURVEY="Maralinga_Far_West_Coast_survey"
SURVEY_STATIONS="['SA176','SA177','SA178','SA178_2','SA179','SA180','SA181','SA182','SA183','SA184','SA26E','SA26W','SA26W_2','SA27M','SA27M_2','SA28','SA29','SA30M','SA31','SA32','SA33','SA34','SA35','SA36','SA37','SA38','SA38_2','SA39','SA40','SA40_2','SA41','SA42','SA43','SA44','SA45','SA46','SA48','SA49','SA50','SA51','SA52','SA53','SA56','SA57_2','SA58','SA59','SA60','SA61','SA61_2','SA62','SA63_2','SA64','SA65','SA67','SA68','SA70','SA71']"

mpirun -np $PBS_NCPUS python3.10 ./mth5_onefileperstation_mt_metadata_v0.2.2_level0.py "$SURVEY" "$SURVEY_STATIONS" > ./pbs_job_logs/$PBS_JOBID.log

date
