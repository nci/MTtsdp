#!/bin/bash

#PBS -N Maralinga_Far_West_Coast_azimuths
#PBS -q normal
#PBS -P <aa0>
#PBS -l walltime=0:10:00
#PBS -l ncpus=16
#PBS -l mem=64GB
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

SURVEY_NAME="Maralinga_Far_West_Coast_survey"

mpirun -np $PBS_NCPUS python3.12 ./02_extract_translated_azimuth.py "$SURVEY_NAME"

date
