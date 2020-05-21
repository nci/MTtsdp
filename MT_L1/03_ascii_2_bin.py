#!/usr/bin/env python3


import glob
import sys
import numpy as np
from mpi4py import MPI
import numpy
import sys
import subprocess

def run(fnm, log_prename):
    comm = MPI.COMM_WORLD
    mpi_size = comm.Get_size()
    rank     = comm.Get_rank()
    pool = open(fnm, 'r').readlines()
    ncmd = len( pool )
    chunk = int(np.ceil(ncmd/mpi_size) )
    i_start, i_end = chunk*rank, chunk*(rank+1)
    if i_end > ncmd:
        i_end = ncmd
    sub_pool = sorted(pool)[i_start:i_end]
    fp = open('%s_%03d.log' % (log_prename, rank), 'w' )
    for cmd in sub_pool:
        cmd = cmd.strip()
        print(cmd, file=fp, flush=True, end='')
        process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        print(output, error, file=fp, flush=True)
    fp.close()

if __name__ == "__main__":
    ### WA
    run('03_workspace/WA.cmd', '03_workspace/ascii2bin_WA_')
    ### SA
    run('03_workspace/SA.cmd', '03_workspace/ascii2bin_SA_')

