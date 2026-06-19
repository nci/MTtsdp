#!/bin/bash

#PBS -N mt_metadata_spreadsheet_to_json
#PBS -q normal
#PBS -P abc
#PBS -l walltime=0:05:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=1GB
#PBS -l storage=gdata/ab12+gdata/up99


module use /g/data/up99/modulefiles
module load NCI-geophys/26.03

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo

SURVEY_SPREADSHEET="NE_SA_survey_mt_metadata_v0.2.2_downsampled.xlsx"
SURVEY="NE_SA_survey"

python3 ./excel_to_json_mt_metadata_v0.2.2_downsampled.py "$SURVEY_SPREADSHEET" "$SURVEY" > ./pbs_job_logs/$PBS_JOBID.log

date
