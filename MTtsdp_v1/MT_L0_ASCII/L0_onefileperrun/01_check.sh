#!/bin/bash

#PBS -N Maralinga_L0_checks
#PBS -q normal
#PBS -P abc
#PBS -l walltime=0:05:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=10GB
#PBS -l storage=scratch/abc+gdata/up99

### current stable version of NCI-geophys module
module use /g/data/up99/modulefiles
module load NCI-geophys/24.08

date

cd ${PBS_O_WORKDIR}

echo
pwd

echo
python3.10 ./01_check.py > 01_workspace/01.log

date
