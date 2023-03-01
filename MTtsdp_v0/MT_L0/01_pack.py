#!/usr/bin/env python3

from mpi4py import MPI
from glob import glob
from os import path
import os
import shutil
import datetime

def check_file(filename, day_name):
    return filename.split('/')[-1].startswith(day_name)

def check_date(day, filename):
    fday = datetime.datetime.strptime(filename.split('/')[-1].split('_')[1][:6], '%y%m%d').timetuple().tm_yday
    return day == str(fday)

def run(mpi_comm, state='SA', state_sub='SA',data_location='00_data_virtual_link/01_data_location/SA/', root_dir= '01_workspace/', log_prefnm = '01_workspace/SA_mpi_log' ):
    ###
    #
    ###
    channels = ['BY', 'BX', 'EX', 'EY','BZ', 'TP', 'ambientTemperature']
    state = state
    state_sub = state_sub
    data_type = 'Packed_Raw_Time_Series_Archive'#'TS'
    period = 'Long_period' # BB or Long_Period
    survey_name = 'SA_AusLAMP_MT_Survey_Musgraves_APY_2016_to_2018'
    frequency = 10
    data_location = data_location
    root_dir = root_dir
    out_location = path.join(root_dir, state, period, survey_name, state_sub, data_type)
    ###
    # 
    ###
    rank = mpi_comm.Get_rank()
    if rank == 0:
        if not os.path.exists(out_location):
            os.makedirs(out_location)
    mpi_comm.Barrier()
    ###
    #
    ###
    exclude = ['Plots', 'config', 'log', 'temp', 'old']
    exclude_folders = ['WA11','WA12','WA28', 'WASA352']
    metadata_dirs = ['config']
    search = path.join(data_location, '*')
    subdirs = sorted([i for i in glob(search) if i.split('/')[-1] not in exclude_folders])
    ###
    # 
    ###
    rank = mpi_comm.Get_rank()
    size = mpi_comm.Get_size()
    chunk = int(len(subdirs) / size) + 1
    i1 = rank*chunk
    i2 = i1+chunk
    logfnm = log_prefnm + '_%03d.txt' % (rank)
    fp = open(logfnm, 'w')
    print('>>> all subdirs: \n\t%s\n\n' % '\n\t'.join(subdirs)  , file=fp, flush=True)
    print('>>> work[%d,%d): \n\t%s\n\n' % (i1, i2, '\n\t'.join(subdirs[i1:i2]) ) , file=fp, flush=True)
    ###
    #
    ###
    for idx, subdir in enumerate(subdirs[i1:i2]):
        if subdir.split('/')[-1].upper() not in exclude_folders:
            site_name = subdir.split('/')[-1].upper()
            output_dir = path.join(out_location, site_name.upper())
            search2 = path.join(subdir, '*')
            vocab_inclusion = subdir + '/' + site_name + '_' + 'RawTimeSeries'
            include = [vocab_inclusion]
            data_dirs = sorted([i for i in glob(search2) if i in include])
            print('>>> (%d/%d) %s' % (idx, i2-i1, subdir) , file=fp, flush=True)
            for folder in data_dirs:
                out_dir_name = folder.split('/')[-1].split('_')[0] #+ '.zip'
                out_dir = path.join(output_dir, out_dir_name)
                print('\t\t%s\n\t\t  ==>%s.zip' % (folder, out_dir), file=fp, flush=True)
                shutil.make_archive(out_dir, 'zip', folder)
    ###
    #
    ###
    mpi_comm.Barrier()

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    ###
    run(comm, state='SA', state_sub='SA',data_location='00_data_virtual_link/01_data_location/SA/', root_dir= '01_workspace/', log_prefnm = '01_workspace/SA_mpi_log')
    ###
    run(comm, state='WA', state_sub='WA',data_location='00_data_virtual_link/01_data_location/WA/', root_dir= '01_workspace/', log_prefnm = '01_workspace/WA_mpi_log')
