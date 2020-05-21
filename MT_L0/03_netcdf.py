#!/usr/bin/env python3

from mpi4py import MPI
from glob import glob
from os import path
import sys
import os
import re
import shutil
import datetime
import itertools
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from datetime import datetime, timedelta
import pytz, datetime
import calendar
from dateutil.parser import parse
import pytz, datetime
import calendar
from dateutil.parser import parse

def check_file(filename, day_name):
    return filename.split('/')[-1].startswith(day_name)

def check_date(day, filename):
    fday = datetime.datetime.strptime(filename.split('/')[-1].split('_')[1][:6], '%y%m%d').timetuple().tm_yday
    return day == str(fday)

def timestamp_to_timeaxis(timestamp,timezone):
    local = pytz.timezone(timezone)
    start_time_local = '20' + timestamp[:2] + '-' + timestamp[2:4] + '-' + timestamp[4:6] + ' ' + timestamp[6:8] + ':' +  timestamp[8:10] + ':' + timestamp[10:12] 
    start_time_obj_local = datetime.datetime.strptime(start_time_local,'%Y-%m-%d %H:%M:%S')
    local_dt = local.localize(start_time_obj_local, is_dst=None)
    utc_dt = local_dt.astimezone(pytz.utc)
    start_time_UTC = utc_dt.strftime('%Y-%m-%d %H:%M:%S')
    time_obj = parse(start_time_UTC)
    unixtime = calendar.timegm(time_obj.timetuple())
    return unixtime

def groupSequence(lst): 
    res = [[lst[0]]] 
    for i in range(1, len(lst)): 
        if lst[i-1]+1 == lst[i]: 
            res[-1].append(lst[i]) 
  
        else: 
            res.append([lst[i]]) 
    return res

def run(mpi_comm, state = 'SA', state_sub = 'SA',
            metadata_dir = '00_data_virtual_link/03_metadata_dir/SA', 
            data_location = '00_data_virtual_link/03_data_location_SA', 
            root_dir = '03_workspace',
            log_prefnm = '03_workspace/SA_mpi_log' ):
    # this is the timezone we are working in
    timezone = "Australia/Adelaide"  
    ### create list of metadata files to add to netCDF files
    metadata_dir = metadata_dir
    exclude_folders = ['WA11','WA12','WA28','WASA352','auslamp_metadata_template.txt']  ### exclude these for now as they are different instruments
    metadata_files = []
    for file in os.listdir(metadata_dir):
            if file.endswith('.txt'):
                if not file.startswith(("WA11","WA12", "WA28", "WASA352", "auslamp_metadata_template.txt")):
                    metadata_files.append(file)
    metadata_files = sorted(metadata_files)        
    metadata_template = 'auslamp_metadata_template.txt'  ### ASEG sample metadata template
    metadata_template_loc = path.join(metadata_dir,metadata_template)
    channels = ['BY', 'BX', 'EX', 'EY','BZ', 'TP', 'ambientTemperature']
    channels_b = ['BY', 'BX', 'EX', 'EY','BZ']
    state = state
    state_sub = state_sub
    data_type = 'Level_0_Concatinated_Time_Series_NetCDF'
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
    exclude = ['Plots', 'config', 'log', 'temp', 'old']
    exclude_folders = ['WA11','WA12','WA28','WASA352']
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
    for idx, subdir in enumerate(subdirs[i1:i2] ):
        if subdir.split('/')[-1].upper() not in exclude_folders:
            print('>>> (%d/%d) %s' % (idx, i2-i1, subdir) , file=fp, flush=True)
            syte_name = subdir.split('/')[-1].upper()
            output_dir = path.join(out_location, syte_name.upper())
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            search2 = path.join(subdir, '*')
            data_dirs = [sorted([i for i in glob(search2)])]
            for folder in data_dirs:
                day_name = [eye.split('/')[-1] for eye in folder]
                for fold,day in zip(folder,day_name):
                    out_dir = path.join(output_dir, day)
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    channel_files_ex = sorted([i for i in os.listdir(fold) if i.endswith(("EX"))])
                    channel_files_ey = sorted([i for i in os.listdir(fold) if i.endswith(("EY"))])
                    channel_files_bx = sorted([i for i in os.listdir(fold) if i.endswith(("BX"))])
                    channel_files_by = sorted([i for i in os.listdir(fold) if i.endswith(("BY"))])
                    channel_files_bz = sorted([i for i in os.listdir(fold) if i.endswith(("BZ"))])
                    names = [channel_files_ex[0].split('.')[0] + '.nc']
                    for name in names:
                        out_name = path.join(out_dir,name)
                    for i,j,k,l,m in zip(channel_files_ex,channel_files_ey,channel_files_bx,channel_files_by,channel_files_bz):
                        for metadata in metadata_files:
                            if metadata.split('.')[0] == syte_name:
                                metadata_loc = path.join(metadata_dir,metadata)
                                fd = open(metadata_loc,'r')
                                lines = [line for line in fd]
                                lines = [line.rstrip('\n') for line in lines]
                                fd2 = open(metadata_template_loc,'r')
                                auslamp_metadata = [line for line in fd2]
                                auslamp_metadata = [line.rstrip('\n') for line in auslamp_metadata]
                                auslamp_metadata_no_equal = [s.strip(' =') for s in auslamp_metadata]
                                lines_and_metadata = []
                                for f,b in zip(auslamp_metadata_no_equal, lines):
                                    if b != "NA":
                                        lines_and_metadata.append([f,b])

                                for f,b in zip([i[0] for i in lines_and_metadata], [i[1] for i in lines_and_metadata]):
                                    globals()[f] = b    
                        EX1 = []
                        EY1 = []
                        BX1 = []
                        BY1 = []
                        BZ1 = []
                        df_bx = open(path.join(fold,k))
                        for line in df_bx:
                            BX1.append(line)        
                        df_by = open(path.join(fold,l))
                        for line in df_by:
                            BY1.append(line)
                        df_bz = open(path.join(fold,m))
                        for line in df_bz:
                            BZ1.append(line)
                        df_ex = open(path.join(fold,i))
                        for line in df_ex:
                            EX1.append(line)   
                        df_ey = open(path.join(fold,j))
                        for line in df_ey:
                            EY1.append(line)

                        timestamp = i.split('_')[1].split('.')[0]
                        unixtime = timestamp_to_timeaxis(timestamp,timezone)
                        t1 = float(unixtime)
                        number_of_samples = len(EX1)
                        seconds = number_of_samples/frequency
                        t2 = t1 + seconds
                        time1 = np.linspace(t1,t2,num = len(EX1))

                        dataset_P = Dataset(out_name, 'w', format='NETCDF4')
                        bx = dataset_P.createDimension('bx', len(BX1))
                        by = dataset_P.createDimension('by', len(BY1))
                        bz = dataset_P.createDimension('bz', len(BZ1))
                        ex = dataset_P.createDimension('ex', len(EX1))
                        ey = dataset_P.createDimension('ey', len(EY1))
                        time = dataset_P.createDimension('time', len(time1))

                        ex = dataset_P.createVariable('ex',np.float32,('ex',))
                        ex.units = 'mV'
                        ex.long_name = 'electric field in N-S orientation'
                        ex.standard_name = 'electric field'
                        try:
                            ex.sampling_rate_Hz = Ex_NS_sample_rate_samples_per_second
                        except Exception:
                            pass    
                        try:
                            ex.NS_electrode_dipole_length_m = Ex_NS_electrode_dipole_length_m
                        except Exception:
                            pass
                        try:    
                            ex.NS_electrode_azimuth_dec_deg = Ex_NS_electrode_azimuth_dec_deg
                        except Exception:
                            pass
                        try:
                            ex.NS_field_gain = Ex_NS_field_gain
                        except Exception:
                            pass
                        try:
                            ex.North_length_m = Ex_North_length_m
                        except Exception:
                            pass
                        try:
                            ex.North_electrode_number = Ex_North_electrode_number
                        except Exception:
                            pass
                        try:
                            ex.North_resistance_to_ground_kOhms = Ex_North_resistance_to_ground_kOhms
                        except Exception:
                            pass
                        try:
                            ex.South_length_m = Ex_South_length_m
                        except Exception:
                            pass
                        try:
                            ex.South_electrode_number = Ex_South_electrode_number
                        except Exception:
                            pass
                        try:
                            ex.South_resistance_to_ground_kOhms = Ex_South_resistance_to_ground_kOhms
                        except Exception:
                            pass
                        try:
                            ex.ground_center_electrode_number = ground_center_electrode_number
                        except Exception:
                            pass
                        try:
                            ex.N_S_center_resistance_to_ground_kOhms = Ex_N_S_center_resistance_to_ground_kOhms
                        except Exception:
                            pass
                        try:
                            ex.NS_AC_mV = Ex_NS_AC_mV
                        except Exception:
                            pass
                        try:
                            ex.NS_DC_mV = Ex_NS_DC_mV
                        except Exception:
                            pass
                        try:
                            ex.MT_recorder_Ex_NS_voltage_mV = MT_recorder_Ex_NS_voltage_mV
                        except Exception:
                            pass
                        ex[:] = EX1

                        ey = dataset_P.createVariable('ey',np.float32,('ey',))
                        ey.units = 'mV'
                        ey.long_name = 'electric field in E-W orientation'
                        ey.standard_name = 'electric field'
                        try:
                            ey.sampling_rate_Hz = Ey_EW_sample_rate_samples_per_second
                        except Exception:
                            pass
                        try:
                            ey.EW_electrode_dipole_length_m = Ey_EW_electrode_dipole_length_m
                        except Exception:
                            pass
                        try:
                            ey.EW_electrode_azimuth_dec_deg = Ey_EW_electrode_azimuth_dec_deg
                        except Exception:
                            pass
                        try:    
                            ey.EW_field_gain = Ey_EW_field_gain
                        except Exception:
                            pass
                        try:    
                            ey.East_length_m = Ey_East_length_m
                        except Exception:
                            pass
                        try:    
                            ey.East_electrode_number = Ey_East_electrode_number
                        except Exception:
                            pass
                        try:    
                            ey.East_resistance_to_ground_kOhms = Ey_East_resistance_to_ground_kOhms
                        except Exception:
                            pass
                        try:    
                            ey.West_length_m = Ey_West_length_m
                        except Exception:
                            pass
                        try:    
                            ey.West_electrode_number = Ey_West_electrode_number
                        except Exception:
                            pass
                        try:    
                            ey.West_resistance_to_ground_kOhms = Ey_West_resistance_to_ground_kOhms
                        except Exception:
                            pass
                        try:    
                            ey.ground_center_electrode_number = ground_center_electrode_number
                        except Exception:
                            pass
                        try:    
                            ey.E_W_center_resistance_to_ground_kOhms = Ey_E_W_center_resistance_to_ground_kOhms
                        except Exception:
                            pass
                        try:    
                            ey.EW_AC_mV = Ey_EW_AC_mV
                        except Exception:
                            pass
                        try:    
                            ey.EW_DC_mV = Ey_EW_DC_mV
                        except Exception:
                            pass
                        try:    
                            ey.MT_recorder_Ey_EW_voltage_mV = MT_recorder_Ey_EW_voltage_mV
                        except Exception:
                            pass
                        ey[:] = EY1

                        bx = dataset_P.createVariable('bx',np.float32,('bx',))
                        bx.units = 'microVolts'
                        bx.long_name = 'instrument recorded magnetic field in X direction (microVolts)'
                        bx.standard_name = 'magnetic field'
                        try:
                            bx.sampling_rate_Hz = Bx_sampling_rate
                        except Exception:
                            pass
                        try:
                            bx.azimuth_degrees = Bx_azimuth_degrees
                        except Exception:
                            pass
                        try:
                            bx.coil_number = Bx_coil_number
                        except Exception:
                            pass
                        bx[:] = BX1

                        by = dataset_P.createVariable('by',np.float32,('by',))
                        by.units = 'microVolts'
                        by.long_name = 'instrument recorded magnetic field in Y direction (microVolts)'
                        by.standard_name = 'magnetic field'
                        try:
                            by.sampling_rate_Hz = By_sampling_rate
                        except Exception:
                            pass
                        try:
                            by.azimuth_degrees = By_azimuth_degrees
                        except Exception:
                            pass
                        try:
                            by.coil_number = By_coil_number
                        except Exception:
                            pass
                        by[:] = BY1

                        bz = dataset_P.createVariable('bz',np.float32,('bz',))
                        bz.units = 'microVolts'
                        bz.long_name = 'instrument recorded magnetic field in Z direction (microVolts)'
                        bz.standard_name = 'magnetic field'
                        try:
                            bz.sampling_rate_Hz = Bz_sampling_rate
                        except Exception:
                            pass
                        try:
                            bz.azimuth_degrees = Bz_azimuth_degrees
                        except Exception:
                            pass
                        try:
                            bz.coil_number = Bz_coil_number
                        except Exception:
                            pass
                        bz[:] = BZ1

                        time = dataset_P.createVariable('time',np.float32,('time',))
                        time.units = 'seconds'
                        time.long_name = 'time'
                        time.standard_name = 'time'
                        try:
                            time.sampling_rate_Hz = Bz_sampling_rate
                        except Exception:
                            pass
                        time[:] = time1


                        try:
                            dataset_P.title = deployment_survey_prospect_name
                        except Exception:
                            pass                     
                        try:
                            dataset_P.summary = 'Long period time series for site'+ ' ' + site_number + ' ' + 'of the' + ' ' + deployment_survey_prospect_name + ' ' + 'survey'
                        except Exception:
                            pass
                        try:
                            dataset_P.source = operator
                        except Exception:
                            pass
                        try:
                            dataset_P.date_created = current_time
                        except Exception:
                            pass
                        try:
                            dataset_P.Conventions = 'ACDD-1.3'
                        except Exception:
                            pass                                    
                            #dataset_P.metadata_link = ' https://geoenetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/fxxxxx'
                        try:
                            dataset_P.license = 'Creative Commons Attribution 4.0 International (CC BY 4.0)'
                        except Exception:
                            pass                                
                            ### Deployment metadata
                        try:
                            dataset_P.site_number = site_number
                        except Exception:
                            pass
                        try:
                            dataset_P.site_name = site_name
                        except Exception:
                            pass
                        try:
                            dataset_P.location = location
                        except Exception:
                            pass
                        try:
                            dataset_P.deployment_survey_prospect_name = deployment_survey_prospect_name
                        except Exception:
                            pass
                        try:
                            dataset_P.institution_agency = institution_agency
                        except Exception:
                            pass
                        try:
                            dataset_P.AUSGRID_number = AUSGRID_number
                        except Exception:
                            pass
                        try:
                            dataset_P.operator = operator
                        except Exception:
                            pass
                        try:
                            dataset_P.recorder_latitude_dec_deg = recorder_latitude_dec_deg
                        except Exception:
                            pass
                        try:
                            dataset_P.recorder_longitude_dec_deg = recorder_longitude_dec_deg
                        except Exception:
                            pass
                        try:
                            dataset_P.recorder_GPS_projection_system = recorder_GPS_projection_system
                        except Exception:
                            pass
                        try:
                            dataset_P.field_GPS_latitude_dec_deg = field_GPS_latitude_dec_deg
                        except Exception:
                            pass
                        try:
                            dataset_P.field_GPS_longitude_dec_deg = field_GPS_longitude_dec_deg
                        except Exception:
                            pass
                        try:
                            dataset_P.GPS_easting_AMG = GPS_easting_AMG
                        except Exception:
                            pass
                        try:
                            dataset_P.GPS_northing_AMG = GPS_northing_AMG
                        except Exception:
                            pass
                        try:
                            dataset_P.AMP_zone = AMG_zone
                        except Exception:
                            pass
                        try:
                            dataset_P.GPS_projection_system = GPS_projection_system
                        except Exception:
                            pass
                        try:
                            dataset_P.altitude_or_elevation_from_GPS_or_instrument = altitude_or_elevation_from_GPS_or_instrument
                        except Exception:
                            pass
                        try:
                            dataset_P.deployment_date = deployment_date
                        except Exception:
                            pass
                        try:
                            dataset_P.LOCAL_start_hour = LOCAL_start_hour
                        except Exception:
                            pass
                        try:
                            dataset_P.LOCAL_start_minute = LOCAL_start_minute
                        except Exception:
                            pass
                        try:
                            dataset_P.LOCAL_start_second = LOCAL_start_second
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_date_from_instrument = UTC_date_from_instrument
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_year = UTC_start_year
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_month = UTC_start_month
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_day = UTC_start_day
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_time = UTC_start_time
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_time_hour = UTC_start_time_hour
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_time_minute = UTC_start_time_minute
                        except Exception:
                            pass
                        try:
                            dataset_P.UTC_start_time_second = UTC_start_time_second
                        except Exception:
                            pass
                        try:
                            dataset_P.deployment_Julian_day = deployment_Julian_day
                        except Exception:
                            pass
                        try:
                            dataset_P.recording_method = recording_method
                        except Exception:
                            pass
                        try:
                            dataset_P.MT_recorder_type_model = MT_recorder_type_model
                        except Exception:
                            pass
                        try:
                            dataset_P.magnetometer_type_model = magnetometer_type_model
                        except Exception:
                            pass
                        try:
                            dataset_P.electrode_type_model = electrode_type_model
                        except Exception:
                            pass
                        try:
                            dataset_P.power_source_type_or_model = power_source_type_or_model
                        except Exception:
                            pass
                        try:
                            dataset_P.data_confidentiality = data_confidentiality
                        except Exception:
                            pass
                        try:
                            dataset_P.north_reference = north_reference
                        except Exception:
                            pass
                        try:
                            dataset_P.drift_calculation_end_time_GPS = drift_calculation_end_time_GPS
                        except Exception:
                            pass
                        try:
                            dataset_P.drift_calculation_end_time_instrument = drift_calculation_end_time_instrument
                        except Exception:
                            pass
                        try:
                            dataset_P.instrument_drift_in_seconds = instrument_drift_in_seconds
                        except Exception:
                            pass
                        try:
                            dataset_P.battery_voltage_no_load_V = battery_voltage_no_load_V
                        except Exception:
                            pass
                        try:
                            dataset_P.battery_voltage_on_load_V = battery_voltage_on_load_V
                        except Exception:
                            pass
                        try:
                            dataset_P.battery_voltage_MT_data_logger_V = battery_voltage_MT_data_logger_V
                        except Exception:
                            pass
                        try:
                            dataset_P.MT_box_case_number = MT_box_case_number
                        except Exception:
                            pass
                        try:
                            dataset_P.MT_recorder_data_logger_number = MT_recorder_data_logger_number
                        except Exception:
                            pass
                        try:
                            dataset_P.MT_recorder_interface_box_number = MT_recorder_interface_box_number
                        except Exception:
                            pass
                        try:
                            dataset_P.data_storage_device_number = data_storage_device_number
                        except Exception:
                            pass
                        try:
                            dataset_P.magnetometer_type = magnetometer_type
                        except Exception:
                            pass
                        try:
                            dataset_P.data_file_length_minutes_per_file = data_file_length_minutes_per_file
                        except Exception:
                            pass
                        try:
                            dataset_P.NS_EW_resistance_kOhms = NS_EW_resistance_kOhms
                        except Exception:
                            pass
                        try:
                            dataset_P.NS_EW_voltage_mV = NS_EW_voltage_mV
                        except Exception:
                            pass
                        try:
                            dataset_P.deployment_comments = deployment_comments
                        except Exception:
                            pass
                        try:
                            dataset_P.photo_taken_deployment = photo_taken_deployment
                        except Exception:
                            pass
                        try:
                            dataset_P.photo_identification_numbers_deployment = photo_identification_numbers_deployment
                        except Exception:
                            pass
                            #dataset_P.deployment_sheet_scan = deployment_sheet_scan
                        try:
                            dataset_P.calculated_date_for_retrieval = calculated_date_for_retrieval
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_date_local = retrieval_date_local
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_local_start_year = retrieval_local_start_year
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_local_start_month = retrieval_local_start_month
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_local_start_day = retrieval_local_start_day
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_time_local = retrieval_time_local
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_local_start_hour = retrieval_local_start_hour
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_local_start_minute = retrieval_local_start_minute
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_local_start_second = retrieval_local_start_second
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_date_UTC = retrieval_date_UTC
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_time_UTC = retrieval_time_UTC
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_Julian_day = retrieval_Julian_day
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_Ex_N_S_Center_resistance_kOhms = retrieval_Ex_N_S_Center_resistance_kOhms
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_Ex_AC_mV = retrieval_Ex_AC_mV
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_Ex_DC_mV = retrieval_Ex_DC_mV
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_Ey_E_W_Center_resistance_kOhms = retrieval_Ey_E_W_Center_resistance_kOhms
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_ey_AC_mV = retrieval_ey_AC_mV
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_Ey_DC_mV = retrieval_Ey_DC_mV
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_NS_EW_resistance_kOhms = retrieval_NS_EW_resistance_kOhms
                        except Exception:
                            pass
                        try:
                            dataset_P.MT_unit_operational = MT_unit_operational
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_battery_voltage_V = retrieval_battery_voltage_V
                        except Exception:
                            pass
                        try:
                            dataset_P.electrode_cable_status = electrode_cable_status
                        except Exception:
                            pass
                        try:
                            dataset_P.mag_cable_status = mag_cable_status
                        except Exception:
                            pass
                        try:
                            dataset_P.retrieval_comments = retrieval_comments
                        except Exception:
                            pass
                        try:
                            dataset_P.days_operated = days_operated
                        except Exception:
                            pass
                        try:
                            dataset_P.photos_taken_retrieval = photos_taken_retrieval
                        except Exception:
                            pass
                        try:
                            dataset_P.photo_identification_numbers_retrieval = photo_identification_numbers_retrieval
                        except Exception:
                            pass
                        try:
                            dataset_P.site_data_quality_retrieval = site_data_quality_retrieval
                        except Exception:
                            pass
                        try:
                            dataset_P.redeployed = redeployed                                
                        except Exception:
                            pass  
                        dataset_P.close()
    mpi_comm.Barrier()

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    ###
    run(comm, state = 'SA', state_sub = 'SA',
            metadata_dir = '00_data_virtual_link/03_metadata_dir/SA', 
            data_location = '00_data_virtual_link/03_data_location_SA', 
            root_dir = '03_workspace',
            log_prefnm = '03_workspace/SA_mpi_log' )
    ###
    run(comm, state = 'WA', state_sub = 'WA',
            metadata_dir = '00_data_virtual_link/03_metadata_dir/WA', 
            data_location = '00_data_virtual_link/03_data_location_WA', 
            root_dir = '03_workspace',
            log_prefnm = '03_workspace/WA_mpi_log' )
