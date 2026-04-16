#!/bin/bash
#PBS -N L0_checks
#PBS -q normal
#PBS -P <aa0>
#PBS -l walltime=0:45:00
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/<xy12>+gdata/up99
#PBS -j oe
#PBS -o /g/data/<xy12>/<abc123>/supplementing_files/pbs_job_logs

export survey_name="Maralinga_Far_West_Coast_survey"

export MT_L0_ASCII="/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/${survey_name}/level_0/Concatenated_Time_Series_ASCII_per_day"
export supplements="/g/data/<xy12>/<abc123>/supplementing_files"


if [ ! -d "${supplements}/${survey_name}" ]; then
    mkdir -p "${supplements}/${survey_name}"
fi

wc -l ${MT_L0_ASCII}/*/*/* | sed  '$d' > ${supplements}/${survey_name}/list_fileperrun.log

### Load NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/25.03

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
python3.12 ./01_check_machine.py 



