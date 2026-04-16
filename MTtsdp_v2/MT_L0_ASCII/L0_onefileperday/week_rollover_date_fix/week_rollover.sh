#!/bin/bash

#PBS -P <aa0>
#PBS -q normalsr
#PBS -l ncpus=104
#PBS -l walltime=00:50:00
#PBS -l mem=500Gb
#PBS -l jobfs=10Gb
#PBS -l storage=gdata/up99+gdata/<xy12>
#PBS -l wd

# load openmpi module
module load openmpi

cd ${PBS_O_WORKDIR}

# execute
mpirun -np $PBS_NCPUS ./week_rollover_date_fix_within_multifile_MPI_t2 SA46-2B SA50-2A SA53-2B

