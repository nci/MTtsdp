#!/bin/bash

#PBS -N mt_metadata_spreadsheet_to_json
#PBS -q normal
#PBS -P abc
#PBS -l walltime=0:05:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=1GB
#PBS -l storage=scratch/abc+gdata/up99


### Need to use NCI-geophys/23.04 for this code

module use /g/data/up99/modulefiles
module load NCI-geophys/23.04

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
python3.9 ./excel_to_json_mt_metadata_v0.2.2.py > ./pbs_job_logs/$PBS_JOBID.log

date
