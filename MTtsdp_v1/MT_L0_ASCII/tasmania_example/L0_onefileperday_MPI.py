from glob import glob
from os import path
import os
import shutil
import itertools
import numpy as np
import time
from mpi4py import MPI

startTime = time.time()

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size



channels = sorted(['BY', 'BX', 'BZ', 'EX', 'EY', 'ambientTemperature'])
state = 'Tasmania'
data_type = 'TS'
level = 'level_0'
phase = 'phase_4'
period = 'Long_Period' # BB or Long_Period
survey_name = 'Mainland_Tasmania_and_Flinders_Island'
frequency = 10

data_location = '/scratch/abc/abc123/AusLAMP_Tasmania/raw_data/phase_4'
outdir = '/scratch/abc/abc123/AusLAMP_Tasmania/AuScope_AusLAMP/'
out_location = path.join(outdir, state, survey_name, level, phase)

exclude = ['ASP', 'Plots', 'config', 'log', 'temp', 'old', 
           'FOUND.000', '$RECYCLE.BIN', 'System\ Volume\ Information/', 
           'System Volume Information', 'merged', 'user_string', 
           'tas124.fig', 'tas124.mat', 'tas124.png','rr221','details.txt',
            'naming_problem.txt']

search = path.join(data_location, '*')
subdirs = sorted([i for i in glob(search) if i.split('/')[-1] not in exclude])

LP_list = ['site_201_LP','site_213_LP','site_221_LP','site_227_LP','site_300_LP','site_315_LP','site_332_LP','site_342_LP']

def write_data_dirs(subdirs, LP_list, exclude, out_location):
    for subdir in subdirs:
        site_name = subdir.split('/')[-1]
        if site_name in LP_list:
            site_name_instrument = 'TAS'+ site_name.split('_')[1]+site_name.split('_')[-1]
        else:
            site_name_instrument = 'TAS'+ site_name.split('_')[-1]
        output_dir = path.join(out_location, site_name_instrument)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)


def search_data_dirs(subdirs, LP_list, exclude, out_location):
    list_of_datadirs = []
    output_dirs = []
    for subdir in subdirs:
        site_name = subdir.split('/')[-1]
        if site_name in LP_list:
            site_name_instrument = 'TAS'+ site_name.split('_')[1]+site_name.split('_')[-1]
        else:
            site_name_instrument = 'TAS'+ site_name.split('_')[-1]
        output_dir = path.join(out_location, site_name_instrument)
        output_dirs.append(output_dir)
        search = path.join(subdir, '*')
        #print(search)
        data_dirs = sorted([i for i in glob(search) if i.split('/')[-1] not in exclude])
        list_of_datadirs.append(data_dirs)
    
    return list_of_datadirs, output_dirs


def create_ascii_files_per_day(combination):
    for directory in combination[0]:
        day_name = directory.split('/')[-1]
        data_files = [i for i in glob(path.join(directory, '*'))]
        out_dir = os.path.join(combination[1], day_name)
        for channel in channels:
            channel_files = [i for i in data_files if i.endswith(channel)]
            if len(channel_files) > 1:
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                channel_files = sorted(channel_files)
                fn = os.path.join(out_dir, channel_files[0].split('/')[-1])
                if not os.path.exists(fn):
                    with open(fn, 'w') as fo:
                        for channel_file in channel_files:
                            with open(channel_file, 'r') as fi:
                                for line in fi:
                                    fo.write(line)


if rank==0:
    if not os.path.exists(out_location):
        os.makedirs(out_location)

    write_data_dirs(subdirs, LP_list, exclude, out_location)


comm.Barrier()

data_dirs, output_dirs = search_data_dirs(subdirs, LP_list, exclude, out_location)

#output_dirs_test = output_dirs[:4]
#data_dirs_test = data_dirs[:4]


for i,combination in enumerate(sorted(itertools.zip_longest(data_dirs,output_dirs))):
    if i%size!=rank:
        continue
    create_ascii_files_per_day(combination)


comm.Barrier()

if rank==0:
    print('The script took {0} seconds !'.format(time.time()-startTime))
