from mpi4py import MPI
import h5py
from mth5.mth5 import MTH5
import numpy as np
from os import path
import os, psutil
import glob
import nc_time_axis
import time

from mt_metadata import timeseries as metadata
from mt_metadata.utils.mttime import MTime
from mt_metadata.timeseries.filters import CoefficientFilter

import json

startTime = time.time()

### define MPI comm, rank and size

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

### define working directories and file paths

work_dir = '/g/data/.../.../workdir'
mth5_output_directory = '/g/data/.../.../outdir'
full_path_to_mth5_files = sorted(glob.glob(mth5_output_directory+"/*")) 
full_path_to_ascii_files = sorted(glob.glob(work_dir+"/*"))
mt_metadata_dir = '/g/data/.../.../mt_metadata_dir'
full_path_to_mt_metadata = sorted(glob.glob(mt_metadata_dir+"/*"))

### define raw time series data channels

raw_data_channels = ['EX','EY','BX','BY','BZ','TP','ambientTemperature']

### define stations to go into mth5 file

stations_all = ['SA225-2','SA227',   'SA242',  'SA243',  'SA245',
                'SA247',  'SA248',   'SA249',  'SA250',  'SA251',  
                'SA252',  'SA26W-2', 'SA270',  'SA271',  'SA272',                
                'SA273',  'SA274-2', 'SA274',  'SA275',  'SA276',
                'SA277',  'SA293-2', 'SA294',  'SA295',  'SA296', 
                'SA297',  'SA298',   'SA300',  'SA301',  'SA319',                            
                'SA320-2', 'SA321',  'SA322',  'SA323', 
                'SA324',  'SA325-2', 'SA325',  'SA326N', 'SA326S',
                'SA344',  'SA344-2', 'SA345',  'SA346',  'SA347',  
                'SA348',  'SA349',   'SA350',  'SA351', 'WA10',            ### 49 single run SA stations
                'WA13',    'WA14',   'WA15',   'WA26', 'WA27',
                'WA29',    'WA30',   'WA31',   'WA42',
                'WA43',   'WA44',    'WA45',   'WA46',   'WA47',
                'WA54',   'WA55',    'WA56',   'WA57',   'WA58',
                'WA60',   'WA61',    'WA62',   'WA63',   'WA64',
                'WA65',   'WA66',    'WA67',   'WA68',   'WA69',
                'WA70',   'WA71',    'WA72',   'WA73',   'WA74',
                'WA75',   'WANT19',  'WANT38', 'WANT45', 'WASA302', 
                'WASA327',                          ### 41 single run WA stations 
                'SA246',  'SA299',   'SA324-2']     ### 3 stations with multiple runs
                                   
                                      
### define survey name and run number for stations with a single run

survey_name = "AusLAMP_Musgraves"
run_number = "001"

### define coefficient filters for Earth Data Logger processing 

gainE = CoefficientFilter(units_in="microvolts", units_out="microvolts", gain=1.0, name="gain_E")
gainB = CoefficientFilter(units_in="microvolts", units_out="microvolts", gain=1.0, name="gain_B")
gainEonly = CoefficientFilter(units_in="microvolts", units_out="microvolts", gain=10.0, name="gain_Eonly")
bz_adjustment = CoefficientFilter(units_in="microvolts",units_out="microvolts",gain=2.2, name="bz_adjustment")

### define functions used in processing

def make_mth5_dir(mth5_output_directory):
### makes mth5 output directory if it doesn't already exist
    try:
        os.makedirs(mth5_output_directory)
    except FileExistsError:
        # directory already exists
        print('directory already exists!')
        pass

def remove_existing_mth5_files(mth5_files):
### removes previously generated mth5_files
    for mth5_file in mth5_files:
        if path.exists(mth5_file):
            os.unlink(mth5_file)
            print("INFO: Removed existing file {}".format(mth5_file))
        else:
            print("File does not exist")


def channel_line_count(channels):
### counts number of lines in ASCII electromagnetic time series files
    EX = [file for file in channels if file.endswith('EX')]
    EY = [file for file in channels if file.endswith('EY')]
    BX = [file for file in channels if file.endswith('BX')]
    BY = [file for file in channels if file.endswith('BY')]
    BZ = [file for file in channels if file.endswith('BZ')]

    count_ex = sum(1 for line in open(EX[0]))
    count_ey = sum(1 for line in open(EY[0]))
    count_bx = sum(1 for line in open(BX[0]))
    count_by = sum(1 for line in open(BY[0]))
    count_bz = sum(1 for line in open(BZ[0]))

    return count_ex, count_ey, count_bx, count_by, count_bz


def channel_data_extraction(channels):
### extracts electromagnetic time series from ASCII generated files
    EX = [file for file in channels if file.endswith('EX')]
    EY = [file for file in channels if file.endswith('EY')]
    BX = [file for file in channels if file.endswith('BX')]
    BY = [file for file in channels if file.endswith('BY')]
    BZ = [file for file in channels if file.endswith('BZ')]
    with open(EX[0], 'r') as file:
        EX1 = file.read().splitlines()
        ex_ts = np.array(EX1).astype(np.int32)
    with open(EY[0], 'r') as file:
        EY1 = file.read().splitlines()
        ey_ts = np.array(EY1).astype(np.int32)
    with open(BX[0], 'r') as file:
        BX1 = file.read().splitlines()
        bx_ts = np.array(BX1).astype(np.int32)
    with open(BY[0], 'r') as file:
        BY1 = file.read().splitlines()
        by_ts = np.array(BY1).astype(np.int32)
    with open(BZ[0], 'r') as file:
        BZ1 = file.read().splitlines()
        bz_ts = np.array(BZ1).astype(np.int32)

    return ex_ts, ey_ts, bx_ts, by_ts, bz_ts


def create_mth5_group_station_run_channel(station,mt_metadata):
### creates Level 0 mth5 file with associated mt_metadata
    with open(mt_metadata[0], 'r') as json_file:
        json_load = json.load(json_file)

    ### extract survey metadata from JSON metadata files
    survey_dict = json_load['survey']
    survey_group.metadata.from_dict(survey_dict)
    survey_group.write_metadata()

    ### add in filters that are used for Level 1 processing 
    survey_group.filters_group.add_filter(gainE)
    survey_group.filters_group.add_filter(gainB)
    survey_group.filters_group.add_filter(gainEonly)
    survey_group.filters_group.add_filter(bz_adjustment)

    ### extract station metadata from JSON metadata files
    station_dict = json_load['station']
    add_station = m.add_station(station, survey=survey_name)
    add_station.metadata.from_dict(station_dict)
    add_station.write_metadata()

    channels = []
    for file in full_path_to_ascii_files:
        if station in file:
            channels.append(file)
        else:
            continue
    if len(channels) == len(raw_data_channels):
    ### for stations with a single run, extract run, ex, ey, bx, by, bz metadata from associated station metadata JSON files
        run_dict = json_load['run']
        ex_dict = json_load['electric_ex']
        ey_dict = json_load['electric_ey']
        bx_dict = json_load['magnetic_bx']
        by_dict = json_load['magnetic_by']
        bz_dict = json_load['magnetic_bz']

        ### add in run and run metadata to mth5 file
        add_run = m.add_run(station, run_number, survey=survey_name) 
        add_run.metadata.from_dict(run_dict)
        add_run.write_metadata()

        ### added in this code as "channels_recorded_electric" wouldn't populate in the mth5 file with the metadata.from_dict function.
        ### Need to ask Jared about this.
        channels_recorded_electric = run_dict['channels_recorded_electric']
        add_run.metadata.channels_recorded_electric = channels_recorded_electric

        ### extract electromagnetic time series from Level 0 ASCII files
        ex_ts,ey_ts,bx_ts,by_ts,bz_ts = channel_data_extraction(channels)
        
        ### add ex,ey,bx,by,bz files and associated metadata to mth5 file
        ex = m.add_channel(station, run_number, "ex", "electric", ex_ts, survey=survey_name)
        ey = m.add_channel(station, run_number, "ey", "electric", ey_ts, survey=survey_name)
        bx = m.add_channel(station, run_number, "bx", "magnetic", bx_ts, survey=survey_name)
        by = m.add_channel(station, run_number, "by", "magnetic", by_ts, survey=survey_name)
        bz = m.add_channel(station, run_number, "bz", "magnetic", bz_ts, survey=survey_name)
        
        ex.metadata.from_dict(ex_dict)
        ex.write_metadata()
        ey.metadata.from_dict(ey_dict)
        ey.write_metadata()
        bx.metadata.from_dict(bx_dict)
        bx.write_metadata()
        by.metadata.from_dict(by_dict)
        by.write_metadata()
        bz.metadata.from_dict(bz_dict)
        bz.write_metadata()
        
        #######################################################################################
        # Note: the resizing of datasets causes h5py parallel to fail when running using      #
        # the mpio driver. A workaround is to create 'zeros' arrays of size count_<xx> (see   #
        # above).                                                                             #
        #                                                                                     #
        # ex = m.add_channel(station, run_number, "ex", "electric", None, survey=survey_name) #
        # ey = m.add_channel(station, run_number, "ey", "electric", None, survey=survey_name) #
        # bx = m.add_channel(station, run_number, "bx", "magnetic", None, survey=survey_name) #
        # by = m.add_channel(station, run_number, "by", "magnetic", None, survey=survey_name) #
        # bz = m.add_channel(station, run_number, "bz", "magnetic", None, survey=survey_name) #
        #                                                                                     #
        # ex.hdf5_dataset.resize((count_ex,))                                                 #
        # ey.hdf5_dataset.resize((count_ey,))                                                 #
        # bx.hdf5_dataset.resize((count_bx,))                                                 #
        # by.hdf5_dataset.resize((count_by,))                                                 #  
        # bz.hdf5_dataset.resize((count_bz,))                                                 #
        #######################################################################################
        
        print(station)
        print(time.time())
   
    elif len(channels) > len(raw_data_channels):
        ### same as above for channels with multiple runs
        sort_files = sorted(channels)
        number_of_channels = len(raw_data_channels)
        split_lists = [sort_files[x:x+number_of_channels] for x in range(0, len(sort_files), number_of_channels)]
        mt_metadata_files = sorted(mt_metadata)
        print(split_lists[0])
        for i, (group,mt_meta) in enumerate(zip(split_lists,mt_metadata_files)):
            mrun_number = i+1
            run = "00%i" % mrun_number
            with open(mt_meta, 'r') as json_file:
                json_load = json.load(json_file)
            
            run_dict = json_load['run']
            ex_dict = json_load['electric_ex']
            ey_dict = json_load['electric_ey']
            bx_dict = json_load['magnetic_bx']
            by_dict = json_load['magnetic_by']
            bz_dict = json_load['magnetic_bz']
            
            add_run = m.add_run(station, run, survey=survey_name)            
            add_run.metadata.from_dict(run_dict)
            add_run.write_metadata()

            ex_ts,ey_ts,bx_ts,by_ts,bz_ts = channel_data_extraction(group)

            ex = m.add_channel(station, run, "ex", "electric", ex_ts, survey=survey_name)
            ey = m.add_channel(station, run, "ey", "electric", ey_ts, survey=survey_name)
            bx = m.add_channel(station, run, "bx", "magnetic", bx_ts, survey=survey_name)
            by = m.add_channel(station, run, "by", "magnetic", by_ts, survey=survey_name)
            bz = m.add_channel(station, run, "bz", "magnetic", bz_ts, survey=survey_name)

            ex.metadata.from_dict(ex_dict)
            ex.write_metadata()
            ey.metadata.from_dict(ey_dict)
            ey.write_metadata()
            bx.metadata.from_dict(bx_dict)
            bx.write_metadata()
            by.metadata.from_dict(by_dict)
            by.write_metadata()
            bz.metadata.from_dict(bz_dict)
            bz.write_metadata()
            
            print(station)
            print(time.time())


    elif len(channels) < len(raw_data_channels):
        print('you are likely missing some channels')
        print(station)

    else:
        print('something has gone wrong')


### create mth5 file   

if rank==0:
    make_mth5_dir(mth5_output_directory)
    remove_existing_mth5_files(full_path_to_mth5_files)

comm.Barrier()

for i,station in enumerate(sorted(stations_all)):
    if i%size!=rank:
        continue
    m = MTH5(file_version='0.2.0',shuffle=None,fletcher32=None,compression="gzip",compression_opts=4)
    hdf5_filename = '{}.h5'.format(station)
    h5_fn = mth5_output_directory+'/'+hdf5_filename 
    m.open_mth5(h5_fn, "w") 
    survey_group = m.add_survey(survey_name)
    mt_metadata_file_name = '{}.json'.format(station)
    mt_metadata = [file for file in full_path_to_mt_metadata if file.endswith(mt_metadata_file_name)]
    create_mth5_group_station_run_channel(station,mt_metadata)
    m.close_mth5()


comm.Barrier()
### print total time to run script

if rank==0:
    print('The script took {0} seconds !'.format(time.time()-startTime))
