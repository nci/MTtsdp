import h5py
from mth5.mth5 import MTH5
from mpi4py import MPI
from scipy.signal import decimate
import json
import sys
import numpy as np
import glob
import os
from os import path
import time

from mt_metadata.timeseries.filters import CoefficientFilter
from mt_metadata import timeseries as metadata
from mt_metadata.utils.mttime import MTime

startTime = time.time()

### define MPI comm, rank and size

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

### define working directories and file paths

work_dir = '/scratch/abc/abc123/.../02_workspace/mth5_L0_outdir'  ### directory with the L0 MTH5 files
downsampled_outdir = '/scratch/abc/abc123/.../02_workspace/mth5_L1_outdir'
full_path_to_mth5_files = sorted(glob.glob(work_dir+"/*")) 
full_path_to_downsampled_outdir = sorted(glob.glob(downsampled_outdir+"/*"))
mt_metadata_downsampled_dir = '/scratch/abc/abc123/.../export_downsampled'
full_path_to_downsampled_mt_metadata = sorted(glob.glob(mt_metadata_downsampled_dir+"/*"))


### define the Level 0 mth5 stations that will be processed (i.e., downsampled/rotated/unit conversion).

stations_single_run = ['SA026E','SA026W_2', 'SA029M',  'SA070',  'SA184']    
                                   

stations_multiple_runs = ['SA057_2']

                                      
### define survey name, run number (for stations with a single run) and resampling rate (in Hz)

survey_name = "AusLAMP_Maralinga"
run_number = "001"
resampling_rate_hz = 1

### define coefficient filters for Earth Data Logger processing

gainE = CoefficientFilter(units_in="microvolts", units_out="microvolts", gain=1.0, name="gain_E")
gainB = CoefficientFilter(units_in="microvolts", units_out="microvolts", gain=1.0, name="gain_B")
gainEonly = CoefficientFilter(units_in="microvolts", units_out="microvolts", gain=10.0, name="gain_Eonly")
bz_adjustment = CoefficientFilter(units_in="microvolts",units_out="microvolts",gain=2.2, name="bz_adjustment")

### define functions used in the generation of Level 1 products

def make_mth5_downsampled_outdir(downsampled_outdir):
    try:
        os.makedirs(downsampled_outdir)
    except FileExistsError:
        # directory already exists
        print('directory already exists!')
        pass

def remove_existing_mth5_files(mth5_files):
    for mth5_file in mth5_files:
        if path.exists(mth5_file):
            os.unlink(mth5_file)
            print("INFO: Removed existing file {}".format(mth5_file))
        else:
            print("File does not exist")


def channel_data_extraction_rotation_downsampling_unit_conversion(mth5_file,station,run,resampling_rate):
    ### extracts time series data and metadata from Level 0 MTH5 files and runs rotation, downsampling and unit conversion routine
    gain_e_string = f'Experiment/Surveys/{survey_name}/Filters/coefficient/gain_e'
    gain_b_string = f'Experiment/Surveys/{survey_name}/Filters/coefficient/gain_b'
    gain_eonly_string = f'Experiment/Surveys/{survey_name}/Filters/coefficient/gain_eonly'
    bz_adjustment_string = f'Experiment/Surveys/{survey_name}/Filters/coefficient/bz_adjustment'
    
    gainE = mth5_file[gain_e_string].attrs['gain']
    gainB = mth5_file[gain_b_string].attrs['gain']
    gainEonly = mth5_file[gain_eonly_string].attrs['gain']
    bz_adjustment = mth5_file[bz_adjustment_string].attrs['gain']
 
    ex_string = f'Experiment/Surveys/{survey_name}/Stations/'+station+'/'+run+'/'+'ex'
    ey_string = f'Experiment/Surveys/{survey_name}/Stations/'+station+'/'+run+'/'+'ey'
    bx_string = f'Experiment/Surveys/{survey_name}/Stations/'+station+'/'+run+'/'+'bx'
    by_string = f'Experiment/Surveys/{survey_name}/Stations/'+station+'/'+run+'/'+'by'
    bz_string = f'Experiment/Surveys/{survey_name}/Stations/'+station+'/'+run+'/'+'bz'
    
    ex = mth5_file[ex_string]
    ey = mth5_file[ey_string]
    bx = mth5_file[bx_string]
    by = mth5_file[by_string]
    bz = mth5_file[bz_string]
   
    translated_azimuth_string = f'Experiment/Surveys/{survey_name}/Stations/' + station + '/' + run + '/' + 'bx'
    translated_azimuth = mth5_file[translated_azimuth_string].attrs['translated_azimuth']
    xlength = mth5_file[ex_string].attrs['dipole_length']
    ylength = mth5_file[ey_string].attrs['dipole_length']
    
    meanBx = np.mean(bx[5000:] ) # hard-coded 5000 from the
    meanBy = np.mean(by[5000:] ) # matlab script provided by Bruce.
    z = complex(meanBx, meanBy)  
    RotAngleBo = np.angle(z)
    RotAngleB = RotAngleBo - translated_azimuth/180.0*np.pi
    new_bx = bx*np.cos(RotAngleB) + by*np.sin(RotAngleB)
    new_by = by*np.cos(RotAngleB) - bx*np.sin(RotAngleB)

    ####
    ## Convert voltages of magnetic field into nanoTesla for bx and by
    ####
    # the values are given in microvolts with the bartington sensors having a
    # maximum output of 10V equivalent to 70000nT, plus taking the gain (see
    # gainB and gainE) into consideration
    # Furthermore the data will be downsampled to 1 second due to sensitivity
    # of the Bartington sensor and to reduce the filesize 

    resample_factor = int(10/resampling_rate_hz)

    #bz_adjustment = 2.2; # to fix large Bz values outside the voltage range (a fix by Jingming Duan of GA)

    intermediate_bx = decimate(new_bx, resample_factor) /10**7 * 70000 * gainB
    intermediate_by = decimate(new_by, resample_factor) /10**7 * 70000 * gainB
    intermediate_bz = decimate(bz, resample_factor) /10**7 * 70000 * gainB * bz_adjustment

    intermediate_ex = decimate(ex, resample_factor)/xlength * gainE / gainEonly
    intermediate_ey = decimate(ey, resample_factor)/ylength * gainE / gainEonly
    
    final_ex = np.array(intermediate_ex,dtype='f4')
    final_ey = np.array(intermediate_ey,dtype='f4')
    final_bx = np.array(intermediate_bx,dtype='f4')
    final_by = np.array(intermediate_by,dtype='f4')
    final_bz = np.array(intermediate_bz,dtype='f4')
    
    return final_ex, final_ey, final_bx, final_by, final_bz


def create_mth5_group_station_run_channel_single_run(mth5_file,station,run,resampling_rate_hz,mt_metadata):
    ### creates a Level 1 mth5 file for stations with a single run with the time series data that has had the downsampling/rotation/unit_conversion routine applied
    with open(mt_metadata[0], 'r') as json_file:
        json_load = json.load(json_file)

    survey_dict = json_load['survey']
    survey_group.metadata.from_dict(survey_dict)
    survey_group.write_metadata()
    survey_group.filters_group.add_filter(gainE)
    survey_group.filters_group.add_filter(gainB)
    survey_group.filters_group.add_filter(gainEonly)
    survey_group.filters_group.add_filter(bz_adjustment)

    station_dict = json_load['station']
    add_station = m.add_station(station, survey=survey_name)
    add_station.metadata.from_dict(station_dict)
    add_station.write_metadata()
    run_dict = json_load['run']
    ex_dict = json_load['electric_ex']
    ey_dict = json_load['electric_ey']
    bx_dict = json_load['magnetic_bx']
    by_dict = json_load['magnetic_by']
    bz_dict = json_load['magnetic_bz']

    add_run = m.add_run(station, run_number, survey=survey_name) 
    add_run.metadata.from_dict(run_dict)
    add_run.write_metadata()
 
    ex_ds, ey_ds, bx_ds, by_ds, bz_ds = channel_data_extraction_rotation_downsampling_unit_conversion(mth5_file,station,run,resampling_rate_hz) 

    ex = m.add_channel(station, run_number, "ex", "electric", ex_ds, survey=survey_name)
    ey = m.add_channel(station, run_number, "ey", "electric", ey_ds, survey=survey_name)
    bx = m.add_channel(station, run_number, "bx", "magnetic", bx_ds, survey=survey_name)
    by = m.add_channel(station, run_number, "by", "magnetic", by_ds, survey=survey_name)
    bz = m.add_channel(station, run_number, "bz", "magnetic", bz_ds, survey=survey_name)
                        
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
        

def create_mth5_group_station_run_channel_multiple_runs(mth5_file,station,resampling_rate_hz,mt_metadata):
    ### creates a Level 1 mth5 file for stations with multiple runs with the time series data that has had the downsampling/rotation/unit_conversion routine applied

    with open(mt_metadata[0], 'r') as json_file:
        json_load = json.load(json_file)

    survey_dict = json_load['survey']
    survey_group.metadata.from_dict(survey_dict)
    survey_group.write_metadata()
    survey_group.filters_group.add_filter(gainE)
    survey_group.filters_group.add_filter(gainB)
    survey_group.filters_group.add_filter(gainEonly)
    survey_group.filters_group.add_filter(bz_adjustment)
    
    station_string = f'Experiment/Surveys/{survey_name}/Stations' + '/' + station
    groups = []
    for group in mth5_file[station_string]:
        if group != 'Transfer_Functions':
            groups.append(group)
    groups = sorted(groups)
    for group,mt_meta in zip(groups,mt_metadata):
            run = "%s" % group
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

            ex_ds, ey_ds, bx_ds, by_ds, bz_ds = channel_data_extraction_rotation_downsampling_unit_conversion(mth5_file,station,run,resampling_rate_hz) 

            ex = m.add_channel(station, run, "ex", "electric", ex_ds, survey=survey_name)
            ey = m.add_channel(station, run, "ey", "electric", ey_ds, survey=survey_name)
            bx = m.add_channel(station, run, "bx", "magnetic", bx_ds, survey=survey_name)
            by = m.add_channel(station, run, "by", "magnetic", by_ds, survey=survey_name)
            bz = m.add_channel(station, run, "bz", "magnetic", bz_ds, survey=survey_name)
                        
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


### create Level 1 mth5 files with associated mt_metadata

if rank==0:
    make_mth5_downsampled_outdir(downsampled_outdir)
    remove_existing_mth5_files(full_path_to_downsampled_outdir)

comm.Barrier()

for i,station in enumerate(sorted(stations_single_run)):
    if i%size!=rank:
        continue
    station_l0_timeseries = [file for file in full_path_to_mth5_files if file.endswith('{}.h5'.format(station))]
    mth5_file_l0 = h5py.File(station_l0_timeseries[0],'r')
    m = MTH5(file_version='0.2.0',shuffle=None,fletcher32=None,compression="gzip",compression_opts=4)
    hdf5_filename = '{}.h5'.format(station)
    h5_fn = downsampled_outdir+'/'+hdf5_filename
    m.open_mth5(h5_fn, "w") 
    survey_group = m.add_survey(survey_name)
    mt_metadata_file_name = '{}.json'.format(station)
    mt_metadata = [file for file in full_path_to_downsampled_mt_metadata if file.endswith(mt_metadata_file_name)]
    create_mth5_group_station_run_channel_single_run(mth5_file_l0,station,run_number,resampling_rate_hz, mt_metadata)
    m.close_mth5()
    mth5_file_l0.close()

for i,station in enumerate(sorted(stations_multiple_runs)):
    if i%size!=rank:
        continue
    station_l0_timeseries = [file for file in full_path_to_mth5_files if file.endswith('{}.h5'.format(station))]
    mth5_file_l0 = h5py.File(station_l0_timeseries[0],'r')
    m = MTH5(file_version='0.2.0',shuffle=None,fletcher32=None,compression="gzip",compression_opts=4)
    hdf5_filename = '{}.h5'.format(station)
    h5_fn = downsampled_outdir+'/'+hdf5_filename
    m.open_mth5(h5_fn, "w")
    survey_group = m.add_survey(survey_name)
    mt_metadata_file_name = '{}.json'.format(station)
    mt_metadata = sorted([file for file in full_path_to_downsampled_mt_metadata if file.endswith(mt_metadata_file_name)])
    with open(mt_metadata[0],'r') as station_json:
        json_load = json.load(station_json)
    station_dict = json_load['station']
    add_station = m.add_station(station, survey=survey_name)
    add_station.metadata.from_dict(station_dict)
    add_station.write_metadata()
    create_mth5_group_station_run_channel_multiple_runs(mth5_file_l0, station, resampling_rate_hz, mt_metadata)
    m.close_mth5()
    mth5_file_l0.close()



comm.Barrier()
### print total time to run script

if rank==0:
    print('The script took {0} seconds !'.format(time.time()-startTime))


