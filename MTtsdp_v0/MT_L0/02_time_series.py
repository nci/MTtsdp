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

def run(mpi_comm, state = 'SA', state_sub = 'SA', data_location='00_data_virtual_link/02_data_location/SA/', root_dir= '02_workspace/', log_prefnm = '02_workspace/SA_mpi_log' ):
    ###
    #
    ###
    state, state_sub = state, state_sub
    channels = ['BY', 'BX', 'EX', 'EY','BZ', 'TP', 'ambientTemperature']
    data_type = 'Level_0_Concatinated_Time_Series_ASCII'
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
    

    exclude = ['Plots', 'config', 'log', 'temp', 'old']
    exclude_folders = ['WA11','WA12','WA28', 'WASA352']
    metadata_dirs = ['config']
    search = path.join(data_location, '*')
    subdirs = sorted([i for i in glob(search) if i.split('/')[-1] not in exclude_folders])
    
    exclude_folders = ['WA11','WA12','WA28', 'WASA352']
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
            print('>>> (%d/%d) %s' % (idx, i2-i1, subdir) , file=fp, flush=True)
            site_name = subdir.split('/')[-1].upper()
            output_dir = path.join(out_location, site_name.upper())
            search2 = path.join(subdir, '*')
            vocab_inclusion = subdir + '/' + site_name + '_' + 'RawTimeSeries'
            include = [vocab_inclusion]
            data_dirs = sorted([i for i in glob(search2) if i in include])
            for folder in data_dirs:
                search3 = path.join(folder,'*')
                day_folders = sorted([i for i in glob(search3) if i.split('/')[-1] not in exclude])
                for data in day_folders:
                    day_name = data.split('/')[-1]
                    data_files = [i for i in glob(path.join(data, '*'))]
                    out_dir = path.join(output_dir, day_name)
                    for channel in channels:
                        channel_files = [i for i in data_files if i.endswith(channel)]
                        if len(channel_files) > 1:
                            if not os.path.exists(out_dir):
                                os.makedirs(out_dir)
                            channel_files = sorted(channel_files)
                            base_name = channel_files[0].split('/')[-1]
                            if site_name + '_' not in base_name:
                                base_name = base_name.replace(site_name, site_name + '_')
                            fn = os.path.join(out_dir, base_name)
                            print('\t\t%s' % (fn) , file=fp, flush=True)
                            if not os.path.exists(fn):
                                with open(fn, 'w') as fo:
                                    for channel_file in channel_files:
                                        with open(channel_file, 'r') as fi:
                                            for line in fi:
                                                fo.write(line)
    mpi_comm.Barrier()

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    ###
    run(comm, state='SA', state_sub='SA',data_location='00_data_virtual_link/02_data_location/SA/', root_dir= '02_workspace/', log_prefnm = '02_workspace/SA_mpi_log')
    ###
    run(comm, state='WA', state_sub='WA',data_location='00_data_virtual_link/02_data_location/WA/', root_dir= '02_workspace/', log_prefnm = '02_workspace/WA_mpi_log')