#!/usr/bin/env python3

import glob
import sys
import numpy as np
from mpi4py import MPI
import numpy
import sys
import os

merged_maralinga = '02_workspace/merged_data_maralinga'

try:
    os.makedirs(merged_maralinga)
except:
    # directory already exists
    pass


def mpi_log_print(out, pre=0, msg='empty message', file= sys.stdout, flush=True):
    if out:
        print('\t'*pre + msg, file= file, flush=flush)
class merge_ascii:
    """
    """
    def __init__(self, fnm, lognm = '', exclude_sta = set() ):
        """
        fnm: Filename of the *.log generated by `wc -l *` that is described in README.
        exclude_sta: stations to be excluded.
        """
        self.comm = MPI.COMM_WORLD
        self.mpi_size = self.comm.Get_size()
        self.rank     = self.comm.Get_rank()
        self.mpi_log = open('%s_%03d.log' % (lognm, self.rank), 'w' )
        #
        self.fnm = fnm
        self.exclude_sta = exclude_sta
        #
        self.vol = dict()
        self.rd()
        self.prepare_mpi_job()

    def rd(self):
        """
        Read data and organize data.
        """
        mpi_log_print(True, 0, '# Input data from %s ...' % self.fnm, file=self.mpi_log, flush=True )
        mpi_log_print(True, 1, 'ignore stations: %s' % (', '.join(list(self.exclude_sta))) , file=self.mpi_log, flush=True )
        for line in open(self.fnm):
            n_sample, data_fnm = line.strip().split()
            n_sample = int(n_sample)
            sta, day, tmp = data_fnm.split('/')[-3:]
            if sta in self.exclude_sta:
                continue # get rid of excluded stations
            dat, com = tmp.split('.')
            if sta not in self.vol:
                self.vol[sta] = dict()
            if com not in self.vol[sta]:
                self.vol[sta][com] = dict()
            self.vol[sta][com][day] = {'size': n_sample, 'nm': dat}
            self.rootdir = '/'+'/'.join( open(self.fnm).readline().split()[1].split('/')[:-3] ) + '/'
    def get_path(self, com, sta, day):
        try:
            path = '/'.join( [self.rootdir, sta, day, self.vol[sta][com][day]['nm']] ) +'.'+com
            return path
        except:
            mpi_log_print(True, 0, 'Err: cannot find path for file %s %s %s' % (com, sta, day), file=sys.stderr, flush=True )
    def prepare_mpi_job(self):
        """
        """
        nsta = len( self.vol.keys() )
        chunk = int(np.ceil(nsta/self.mpi_size) )
        i_start, i_end = chunk*self.rank, chunk*(self.rank+1)
        if i_end > nsta:
            i_end = nsta
        sta_lst = sorted(list(self.vol.keys() ) )
        self.proc_sta = sta_lst[i_start: i_end]
        mpi_log_print(True, 0, '# Init mpi job. This node takes care of stations: [%d->%d) (%d/%d)' % (i_start, i_end, len(self.proc_sta ), nsta) , file=self.mpi_log, flush=False)
        mpi_log_print(True, 1, '%s' % ', '.join(self.proc_sta) , file=self.mpi_log, flush=True)
    def merge(self, out_prename, cut_head=1, cut_end=1):
        """
        Merge station by station
        cut_head, cut_end (1, 1) : how many days to remove when merging
        """
        mpi_log_print(True, 0, '# Merging, get rid of starting and end days (%d, %d)...' % (cut_head, cut_end), file=self.mpi_log, flush=True)
        nsta = len(self.proc_sta)
        for ista, sta in enumerate(self.proc_sta):
            for com in self.vol[sta]:
                day_lst = sorted(list(self.vol[sta][com].keys() ) )
                n_days = len(day_lst)
                outfnm = '%s/%s-%s_%s.%s' % (out_prename, sta, day_lst[cut_head], day_lst[cut_end*-1-1], com)
                mpi_log_print(True, 1, '(%d/%d) %s %s ==> %s' % (ista+1, nsta, sta, com, outfnm), file=self.mpi_log, flush=True )
                fp = open(outfnm, 'w')
                for a_day in day_lst[cut_head: n_days-cut_end]:
                    path = self.get_path(com, sta, a_day)
                    mpi_log_print(True, 2, '%s' % path, file=self.mpi_log, flush=True )
                    vol_in = open(path, 'r').readlines()
                    print(''.join(vol_in), end='', file= fp)
                fp.close()



if __name__ == "__main__":
    Maralinga_sta_exclude_set = set( (
        'SA026E','SA026W_2','SA029M','SA070','SA184'   # these are the good stations
    ) )

    ######### Maralinga
    app = merge_ascii('01_workspace/Maralinga.log', '02_workspace/maralinga_merge_', Maralinga_sta_exclude_set)
    app.merge('02_workspace/merged_data_maralinga', 26, 1)

    

