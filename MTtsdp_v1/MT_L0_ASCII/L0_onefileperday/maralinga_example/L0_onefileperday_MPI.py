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
state = 'SA'
data_type = 'TS'
level = 'level_0'
period = 'Long_Period' # BB or Long_Period
survey_name = 'Maralinga_Far_West_Coast_survey'
frequency = 10

data_location = '/scratch/abc/AusLAMP_Maralinga/Maralinga_Far_West_Coast_survey' ### change this to the location of your Earth Data Logger data

outdir = '/scratch/.....' ### change this to the location where you want the new data written to

out_location = path.join(outdir, state, survey_name, level)

#exclude = ['ASP', 'Plots', 'config', 'log', 'temp', 'old', 
#           'FOUND.000', '$RECYCLE.BIN', 'System\ Volume\ Information/', 
#           'System Volume Information', 'merged', 'user_string', 
#           'tas124.fig', 'tas124.mat', 'tas124.png','rr221','details.txt',
#            'naming_problem.txt']

exclude = ['config', 'log', 'temp']


search = path.join(data_location, '*')
subdirs = sorted([i for i in glob(search) if i.split('/')[-1] not in exclude])

def write_data_dirs(subdirs, exclude, out_location):
    for subdir in subdirs:
        site_name = subdir.split('/')[-1]
        output_dir = path.join(out_location, site_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)


def search_data_dirs(subdirs, exclude, out_location):
    list_of_datadirs = []
    output_dirs = []
    for subdir in subdirs:
        site_name = subdir.split('/')[-1]
        output_dir = path.join(out_location, site_name)
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

    write_data_dirs(subdirs, exclude, out_location)


comm.Barrier()

data_dirs, output_dirs = search_data_dirs(subdirs, exclude, out_location)


for i,combination in enumerate(sorted(itertools.zip_longest(data_dirs,output_dirs))):
    if i%size!=rank:
        continue
    create_ascii_files_per_day(combination)


comm.Barrier()

if rank==0:
    print('The script took {0} seconds !'.format(time.time()-startTime))
