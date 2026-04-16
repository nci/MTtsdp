#!/bin/bash
#PBS -N NF_L0_checks
#PBS -q normal
#PBS -P eq93
#PBS -l walltime=1:00:00
#PBS -l ncpus=12
#PBS -l mem=48GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/xr05+gdata/up99
#PBS -j oe
#PBS -o /g/data/xr05/supplementing_files/pbs_job_logs

# export survey_name="Maralinga_Far_West_Coast_survey"
# export survey_name="Eyre_Kangaroo_Island_survey"
# export survey_name="Flinders_Ranges_survey"
# export survey_name="Gawler_survey"
# export survey_name="NE_SA_survey"
export survey_name="North_Flinders_survey"
# export survey_name="SE_SA_survey"


export MT_L0_ASCII="/g/data/xr05/processed_data/South_Australia/${survey_name}/level_0/Concatenated_Time_Series_ASCII_per_day"
export supplements="/g/data/xr05/supplementing_files"


if [ ! -d "${supplements}/${survey_name}" ]; then
    mkdir -p "${supplements}/${survey_name}"
fi

wc -l ${MT_L0_ASCII}/*/*/* | sed  '$d' > ${supplements}/${survey_name}/list_fileperrun.log

### current stable version of NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
python3.10 ./01_check_machine.py 



