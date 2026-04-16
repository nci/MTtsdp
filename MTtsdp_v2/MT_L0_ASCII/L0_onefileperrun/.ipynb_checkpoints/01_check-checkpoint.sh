#!/bin/bash

#PBS -N Maralinga_L0_checks
#PBS -q normal
#PBS -P eq93
#PBS -l walltime=0:10:00
#PBS -l ncpus=6
#PBS -l mem=24GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/xr05+gdata/up99
#PBS -j oe
#PBS -o /g/data/xr05/supplementing_files/pbs_job_logs

MT_L0_ASCII="/g/data/xr05/my80_test_full/outdata/South_Australia/Maralinga_Far_West_Coast_survey/level_0/Concatenated_Time_Series_ASCII_per_day"
supplements="/g/data/xr05/my80_test_full/Maralinga_sample/supplementing_files"
survey="Maralinga"

wc -l ${MT_L0_ASCII}/*/*/* | sed  '$d' > ${supplements}/${survey}_list_fileperrun.log


### current stable version of NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
python3.10 ./01_check.py > ${supplements}/${survey}_filecheck.log

date
