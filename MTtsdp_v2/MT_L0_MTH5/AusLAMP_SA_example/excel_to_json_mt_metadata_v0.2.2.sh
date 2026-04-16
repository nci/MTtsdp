#!/bin/bash

#PBS -N mt_metadata_spreadsheet_to_json
#PBS -q normal
#PBS -P <aa0>
#PBS -l walltime=0:05:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=1GB
#PBS -l storage=gdata/<xy12>+gdata/up99

### Load NCI-geophys/23.04 for this code

module use /g/data/up99/modulefiles
module load NCI-geophys/23.04

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo

SURVEY_SPREADSHEET="Maralinga_mt_metadata_v0.2.2.xlsx"
SURVEY="Maralinga_Far_West_Coast_survey"

python3.9 ./excel_to_json_mt_metadata_v0.2.2.py "$SURVEY_SPREADSHEET" "$SURVEY" > ./pbs_job_logs/$PBS_JOBID.log

date
