import h5py
import shutil
from os import path
import os
import glob
import datetime as dt
import time
from mpi4py import MPI
from mtpy.processing import birrp

startTime = time.time()

### define MPI comm, rank and size

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size


### define working directories, file paths and birrp executable location

frequency = 1  #hertz
bp_workdir = '/g/data/.../.../workdir'
bp_outdir = '/g/data/.../.../outdir'
birrp_location = '/g/data/.../.../birrp'
final_edi_dir = '/g/data/.../.../final_EDI_files'


if rank==0:
    if path.exists(final_edi_dir):
        print('final EDI directory already exists')
    else:
        os.mkdir(final_edi_dir)
        print('created final EDI directory')

### define station combinations. Needs to be in the format [('local_site','remote_site','local_run_number','remote_run_number')
### For example: [('SA250','SA275','001','001'),('SA252','SA275','001','001')] 

station_combinations = [('SA225-2','SA294','001','001'),
                        ('SA227','SA243','001','001'),
                        ('SA242','SA242','001','001'),  ### SA242 uses WAS352 as RR - we don't have this station  
                        ('SA243','SA294','001','001'),
                        ('SA245','SA245','001','001'),
                        ('SA246','SA294','001','001'),  ### SA246 has multiple runs
                        ('SA247','SA294','001','001'),
                        ('SA248','SA275','001','001'),
                        ('SA249','SA275','001','001'),
                        ('SA250','SA275','001','001'),
                        ('SA251','SA275','001','001'),
                        ('SA252','SA275','001','001'), 
                        ('SA26W-2','SA275','001','001'),
                        ('SA270','SA294','001','001'),
                        ('SA271','SA294','001','001'),  ### original script uses SA294 for local and remote?
                        ('SA272','SA294','001','001'),
                        ('SA273','SA275','001','001'),  ### original uses SA273 ex,ey,bz SA251 bx,by SA275 RR
                        ('SA274-2','SA294','001','001'), ### SA274 not used
                        ('SA275','SA252','001','001'),
                        ('SA276','SA275','001','001'),
                        ('SA277','SA275','001','001'),
                        ('SA293-2','SA294','001','001'),
                        ('SA294','SA295','001','001'),
                        ('SA295','SA294','001','001'),
                        ('SA296','SA294','001','001'),
                        ('SA297','SA294','001','001'),
                        ('SA298','SA294','001','001'),
                        ('SA299','SA275','004','001'),  ### SA299 has multiple runs
                        ('SA300','SA275','001','001'),
                        ('SA301','SA275','001','001'),  ### SA271 doesn't overlap - tried SA275 and it worked
                        ('SA319','SA294','001','001'),
                        ('SA320-2','SA320-2','001','001'), ### SA320-2 and SA294 don't have enough overlap
                        ('SA321','SA294','001','001'),
                        ('SA322','SA275','001','001'),
                        ('SA323','SA275','001','001'),
                        ('SA324-2','SA294','001','001'), ### SA324-2 has multiple runs
                        ('SA325-2','SA294','001','001'),
                        ('SA326N','SA275','001','001'),
                        ('SA326S','SA275','001','001'),
                        #('SA344-2','SA344-2','001','001'), ### length of time series section not long enough
                        ('SA345','SA295','001','001'),
                        ('SA346','SA294','001','001'),
                        ('SA347','SA275','001','001'),
                        ('SA348','SA275','001','001'),
                        ('SA349','SA275','001','001'),
                        ('SA350','SA275','001','001'),
                        ('SA351','SA275','001','001'),
                        ('WA10','WA30','001','001'),
                        ('WA13','WA30','001','001'),
                        ('WA14','WA30','001','001'),
                        ('WA15','WA30','001','001'),
                        ('WA26','WA30','001','001'),
                        ('WA27','WA30','001','001'),
                        ('WA29','WA30','001','001'),
                        ('WA30','WA45','001','001'),
                        ('WA31','WA30','001','001'),
                        ('WA42','WA30','001','001'),
                        ('WA43','WA43','001','001'),  ### original uses SA075 as both local and remote??
                        ('WA44','WA30','001','001'),
                        ('WA45','WA46','001','001'),
                        ('WA46','WA45','001','001'),
                        ('WA47','WA30','001','001'),
                        ('WA54','WA61','001','001'),
                        ('WA55','WA61','001','001'),
                        ('WA56','WA61','001','001'),  ### original uses WA56 and WA30 which don't overlap
                        ('WA57','WA58','001','001'),
                        ('WA58','WA57','001','001'),
                        ('WA60','WA61','001','001'),
                        ('WA61','WA71','001','001'),
                        ('WA62','WA57','001','001'),
                        ('WA63','WA57','001','001'),
                        ('WA64','WA57','001','001'),
                        ('WA65','WA57','001','001'),
                        ('WA66','WA61','001','001'),
                        ('WA67','WA61','001','001'),
                        ('WA68','WA57','001','001'),
                        ('WA69','WA61','001','001'),
                        ('WA70','WA61','001','001'),  ### original uses WA70 ex,ey,bz, WA69 bx,by, WA61 bx,by
                        ('WA71','WA61','001','001'),
                        ('WA72','WA61','001','001'),
                        ('WA73','WA57','001','001'),
                        ('WA74','WA57','001','001'),
                        ('WA75','WA57','001','001'),
                        ('WANT19','WA30','001','001'),
                        ('WANT38','WANT38','001','001'),
                        ('WANT45','WA30','001','001'),
                        ('WASA302','WA30','001','001'),
                        ('WASA327','WA30','001','001')] 


def birrp_processing(station_combinations,birrp_workdir,birrp_outdir):
    for i,combination in enumerate(station_combinations):
        if i%size!=rank:
            continue
        station1 = combination[0]
        station2 = combination[1]
        run_station1 = combination[2]
        run_station2 = combination[3]

        file1_orig = birrp_workdir + '/'+ station1 + '.h5'
        file2_orig = birrp_workdir + '/'+ station2 + '.h5'

        os.chdir(birrp_outdir)

        new_directory = birrp_outdir+'/'+station1+'_RR_'+station2

        ### create new_directory for BIRRP processing
        if path.exists(new_directory):
            print('folder already exists')
        else:
            os.mkdir(new_directory)
            print('created new folder')

        os.chdir(new_directory)

        ### copy local and remote reference Level 1 mth5 files into new directory
        shutil.copy(file1_orig,new_directory)
        shutil.copy(file2_orig,new_directory)

        file1 = new_directory + '/'+ station1 + '.h5'
        file2 = new_directory + '/'+ station2 + '.h5'        
        
        ### open mth5 files in read mode using h5py
        f1 = h5py.File(file1,'r')
        f2 = h5py.File(file2,'r')

        ### define survey, station, run, filter strings that will be used for metadata extraction from our L1 mth5 files
        survey_string = 'Experiment/Surveys/AusLAMP_Musgraves'
        
        station1_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station1
        station2_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station2

        f1_run_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station1 + '/' + run_station1 
        f2_run_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station2 + '/' + run_station2 

        ex_run_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station1 + '/' + run_station1 + '/'+'ex'
        ey_run_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station1 + '/' + run_station1 + '/'+'ey'
        
        bx_run_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station1 + '/' + run_station1 + '/'+'bx'
        by_run_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+ station1 + '/' + run_station1 + '/'+'by'

        gain_e_string = 'Experiment/Surveys/AusLAMP_Musgraves/Filters/coefficient/gain_e'
        gain_b_string = 'Experiment/Surveys/AusLAMP_Musgraves/Filters/coefficient/gain_b'
        gain_eonly_string = 'Experiment/Surveys/AusLAMP_Musgraves/Filters/coefficient/gain_eonly'
        bz_adjustment_string = 'Experiment/Surveys/AusLAMP_Musgraves/Filters/coefficient/bz_adjustment'

        ### define run start and end times for local and remote reference stations
        f1_run_starttime = f1[f1_run_string].attrs['time_period.start']
        f1_run_endtime = f1[f1_run_string].attrs['time_period.end']

        f1_run_start = f1_run_starttime[:-6]
        f1_run_end = f1_run_endtime[:-6]

        f2_run_starttime = f2[f2_run_string].attrs['time_period.start']
        f2_run_endtime = f2[f2_run_string].attrs['time_period.end']

        f2_run_start = f2_run_starttime[:-6]
        f2_run_end = f2_run_endtime[:-6]


        ### adapting time stamps to make them useable for later routines 
        if "." in f1_run_start:
            f1_start_element = dt.datetime.strptime(f1_run_start, "%Y-%m-%dT%H:%M:%S.%f")
        else:
            f1_start_element = dt.datetime.strptime(f1_run_start, "%Y-%m-%dT%H:%M:%S")

        if "." in f1_run_end:
            f1_end_element = dt.datetime.strptime(f1_run_end, "%Y-%m-%dT%H:%M:%S.%f")
        else:
            f1_end_element = dt.datetime.strptime(f1_run_end, "%Y-%m-%dT%H:%M:%S")

        if "." in f2_run_start:
            f2_start_element = dt.datetime.strptime(f2_run_start, "%Y-%m-%dT%H:%M:%S.%f")
        else:
            f2_start_element = dt.datetime.strptime(f2_run_start, "%Y-%m-%dT%H:%M:%S")

        if "." in f2_run_end:
            f2_end_element = dt.datetime.strptime(f2_run_end, "%Y-%m-%dT%H:%M:%S.%f")
        else:
            f2_end_element = dt.datetime.strptime(f2_run_end, "%Y-%m-%dT%H:%M:%S")
            
        f1_starttime = dt.datetime.timestamp(f1_start_element)
        f1_endtime = dt.datetime.timestamp(f1_end_element)
        f2_starttime = dt.datetime.timestamp(f2_start_element)
        f2_endtime = dt.datetime.timestamp(f2_end_element)
        
        ### defining the number of samples to skip at the start of the timeseries
        if f1_starttime > f2_starttime:
            number_of_samples_to_skip_start = round(f1_starttime - f2_starttime)
            station_to_skip_start = station2
        elif f1_starttime < f2_starttime:
            number_of_samples_to_skip_start = round(f2_starttime - f1_starttime)
            station_to_skip_start = station1
        elif f1_starttime == f2_starttime:
            number_of_samples_to_skip_start = 0
            station_to_skip_start = None

        if station_to_skip_start == station1:
            number_of_samples_to_skip_station1 = number_of_samples_to_skip_start
        else:
            number_of_samples_to_skip_station1 = 0
    
        if station_to_skip_start == station2:
            number_of_samples_to_skip_station2 = number_of_samples_to_skip_start
        else:
            number_of_samples_to_skip_station2 = 0
    
        if station_to_skip_start == None:
            number_of_samples_to_skip_station1 = 0
            number_of_samples_to_skip_station2 = 0
            
        ### define number of samples to skip at the end of the timeseries and the length_of_timeseries to be used as a BIRRP inupt
        if f1_endtime > f2_endtime:
            number_of_samples_to_skip_end = round(f1_endtime - f2_endtime)
            length_of_timeseries = round(f2_endtime-max(f1_starttime,f2_starttime)) * frequency 
        elif f1_endtime < f2_endtime:
            number_of_samples_to_skip_end = round(f2_endtime - f1_endtime)
            length_of_timeseries = round(f1_endtime-max(f1_starttime,f2_starttime)) * frequency
        elif f1_endtime == f2_endtime:
            number_of_samples_to_skip_end = 0
            length_of_timeseries = round(f1_endtime-max(f1_starttime,f2_starttime)) * frequency
            
        ### extract ex,ey azimuth metadata from mth5 files
        ex_azimuth = f1[ex_run_string].attrs['measurement_azimuth']
        ey_azimuth = f1[ey_run_string].attrs['measurement_azimuth']

        ### define mth5 strings that will be used by BIRRP to extract the L1 time series data 
        ex_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station1+'/'+run_station1+'/'+'ex'+'@'+station1+'.h5'
        ey_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station1+'/'+run_station1+'/'+'ey'+'@'+station1+'.h5'
        bx_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station1+'/'+run_station1+'/'+'bx'+'@'+station1+'.h5'
        by_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station1+'/'+run_station1+'/'+'by'+'@'+station1+'.h5'
        bz_string = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station1+'/'+run_station1+'/'+'bz'+'@'+station1+'.h5'
        bx_string_RR = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station2+'/'+run_station2+'/'+'bx'+'@'+station2+'.h5'
        by_string_RR = 'Experiment/Surveys/AusLAMP_Musgraves/Stations/'+station2+'/'+run_station2+'/'+'by'+'@'+station2+'.h5'
        
        ### Define BIRRP processing parameters
        input_level = '0'
        number_of_output_time_series = '3'
        number_of_input_time_series = '2'
        tbw_for_prolate_data_window = '2'
        data_sample_interval = '1'
        initial_section_length_and_maximum_number_of_sections = '65536,13'
        are_these_values_acceptable = 'y'
        robustness_and_leverage_parameters = '0,0'
        second_stage_coherence_threshold = '0'
        third_varibale_threshold_mode = '0'
        second_stage_coherence_threshold = '0'
        output_filename_root = '{}_RR_{}'.format(station1,station2)
        output_level = '0'
        number_of_data_pieces = '1'
        ar_filter_length = '1'
        file_mode = '0'
        input_mode = '0'
        number_of_points_to_read = '{}'.format(length_of_timeseries)
        number_of_filter_parameters1 = '0'
        data_filename1_ex = ex_string
        number_of_points_to_skip_1 = '{}'.format(number_of_samples_to_skip_station1)
        number_of_filter_parameters2 = '0'
        data_filename1_ey = ey_string
        number_of_points_to_skip_2 = '{}'.format(number_of_samples_to_skip_station1)
        number_of_filter_parameters3 = '0'
        data_filename1_bz = bz_string
        number_of_points_to_skip_3 = '{}'.format(number_of_samples_to_skip_station1)
        number_of_filter_parameters4 = '0'
        data_filename1_bx = bx_string
        number_of_points_to_skip_4 = '{}'.format(number_of_samples_to_skip_station1)
        number_of_filter_parameters5 = '0'
        data_filename1_by = by_string
        number_of_points_to_skip_5 = '{}'.format(number_of_samples_to_skip_station1)
        number_of_filter_parameters6 = '0'
        data_filename2_bx = bx_string_RR
        number_of_points_to_skip_6 =  '{}'.format(number_of_samples_to_skip_station2)
        number_of_filter_parameters7 = '0'
        data_filename2_by = by_string_RR
        number_of_points_to_skip_7 =  '{}'.format(number_of_samples_to_skip_station2)
        rotangle1 = '{},{},180'.format(round(ex_azimuth),round(ey_azimuth))
        rotangle2 = '0,90,0'
        rotangle3 = '0,90,0'


        ### define extra BIRRP parameters required to create EDI files
        max_window_length = '65536'
        n_bisections = '13'
        phi = '180'
        stationlist = [station1]
        toplevelstation = '[%s]' % ', '.join(map(str, stationlist))

        ### define BIRRP string that will be used as input for the BIRRP code
        
        birrp_string = [input_level, 
                number_of_output_time_series,
                number_of_input_time_series, 
                tbw_for_prolate_data_window, 
                data_sample_interval, 
                initial_section_length_and_maximum_number_of_sections, 
                are_these_values_acceptable, 
                robustness_and_leverage_parameters, 
                second_stage_coherence_threshold,
                third_varibale_threshold_mode, 
                second_stage_coherence_threshold, 
                output_filename_root, 
                output_level, 
                number_of_data_pieces,
                ar_filter_length, 
                file_mode, 
                input_mode, 
                number_of_points_to_read,
                number_of_filter_parameters1,
                data_filename1_ex,
                number_of_points_to_skip_1,
                number_of_filter_parameters2,
                data_filename1_ey,
                number_of_points_to_skip_2,
                number_of_filter_parameters3,
                data_filename1_bz,
                number_of_points_to_skip_3,
                number_of_filter_parameters4,
                data_filename1_bx,
                number_of_points_to_skip_4,
                number_of_filter_parameters5,
                data_filename1_by,
                number_of_points_to_skip_5,
                number_of_filter_parameters6,
                data_filename2_bx,
                number_of_points_to_skip_6,
                number_of_filter_parameters7,
                data_filename2_by,
                number_of_points_to_skip_7,
                rotangle1,
                rotangle2,
                rotangle3]
        
        ### write BIRRP script file
        
        birrp_script_name = new_directory+'/'+station1+'_RR_'+station2+'.script'
        
        with open(birrp_script_name,'w') as fp:
            for item in birrp_string:
                fp.write("%s\n" % item)
            print('The birrp script has been written')
            

        ### define and write BIRRP config file that is required to generate an EDI file when using MTpy's J2Edi class ( https://github.com/MTgeophysics/mtpy/blob/217eb5050f5fa161a151be18ed09be105683d0fe/mtpy/processing/birrp.py )
        
        birrp_config_file_name = new_directory+'/'+station1+'_RR_'+station2+'_'+'birrpconfig.cfg'

        birrp_config_string = [toplevelstation,
        'ainuin = ' '{}'.format(robustness_and_leverage_parameters),
        'coherence_threshold = ' '{}'.format(second_stage_coherence_threshold),
        'ilev = '  '{}'.format(input_level),
        'imode = ' '{}'.format(input_mode),
        'jmode = ' '{}'.format(file_mode),
        'max_window_length = ' '{}'.format(max_window_length),
        'n_bisections = ' '{}'.format(n_bisections),
        'n_output_channels = ' '{}'.format(number_of_output_time_series),
        'n_samples = ' '{}'.format(number_of_points_to_read),
        'nar = ' '{}'.format(ar_filter_length),
        'nfil = ' '{}'.format(number_of_filter_parameters1),
        'ninp = ' '{}'.format(number_of_input_time_series),
        'nlev = ' '{}'.format(output_level),
        'nout = ' '{}'.format(number_of_output_time_series),
        'npcs = ' '{}'.format(number_of_data_pieces),
        'nskip = ' '{}'.format(number_of_points_to_skip_7),
        'phi = ' '{}'.format(phi),
        'rr_station = ' '{}'.format(station2),
        'sampling_rate = ' '{}'.format(frequency),
        'station = ' '{}'.format(station1),
        'tbw = ' '{}'.format(tbw_for_prolate_data_window),
        'theta1 = ' '{}'.format(ex_azimuth),
        'theta2 = ' '{}'.format(ey_azimuth)]
             
        with open(birrp_config_file_name,'w') as fp:
            for item in birrp_config_string:
                fp.write("%s\n" % item)
            print('The birrp config file has been written')
        
        ### define survey config file parameters that are required to generate an EDI file when using MTpy's J2Edi class ( https://github.com/MTgeophysics/mtpy/blob/217eb5050f5fa161a151be18ed09be105683d0fe/mtpy/processing/birrp.py )
        
        sampling = f1[f1_run_string].attrs['sample_rate']
        acquired_by =  f1[f1_run_string].attrs['acquired_by.author']
        b_instrument_manufacturer = f1[bx_run_string].attrs['sensor.manufacturer']
        b_instrument_model = f1[bx_run_string].attrs['sensor.model']
        b_logger_type = f1[f1_run_string].attrs['data_logger.manufacturer']
        b_instrument_type = f1[bx_run_string].attrs['sensor.type']

        b_logger_gain = f1[gain_b_string].attrs['gain']
        b_instrument_amplification = f1[bz_adjustment_string].attrs['gain']

        b_xaxis_azimuth = f1[bx_run_string].attrs['measurement_azimuth'] 
        b_yaxis_azimuth = f1[by_run_string].attrs['measurement_azimuth']

        e_logger_type = f1[f1_run_string].attrs['data_logger.type']
        e_logger_gain = f1[gain_e_string].attrs['gain']
        e_instrument_type = f1[f1_run_string].attrs['data_logger.manufacturer']
        e_instrument_amplification = f1[gain_eonly_string].attrs['gain']
        e_xaxis_azimuth = f1[ex_run_string].attrs['measurement_azimuth']
        e_xaxis_length = f1[ex_run_string].attrs['dipole_length']
        e_yaxis_azimuth = f1[ey_run_string].attrs['measurement_azimuth']
        e_yaxis_length = f1[ey_run_string].attrs['dipole_length']

        data_logger_manufacturer = f1[f1_run_string].attrs['data_logger.manufacturer']
        data_logger_type = f1[f1_run_string].attrs['data_logger.type']
        lat = f1[station1_string].attrs['location.latitude']
        lon = f1[station1_string].attrs['location.longitude']
        elevation = f1[station1_string].attrs['location.elevation']
        location = f1[survey_string].attrs['name']
        station = station1
        station_type = 'mt'

        date = f1[survey_string].attrs['time_period.start_date']
        network = f1[survey_string].attrs['project']
        hx = 'N/A'
        hy = 'N/A'
        hz = 'N/A'
        rr_box = 'N/A'
        rr_date = f2[survey_string].attrs['time_period.start_date']
        rr_hx = 'N/A'
        rr_hy = 'N/A'
        rr_lat = f2[station2_string].attrs['location.latitude']
        rr_lon = f2[station2_string].attrs['location.longitude']
        rr_station = station2
             

        ### Define survey config string and write survey config file

        survey_config_file_name = new_directory+'/'+station1+'_RR_'+station2+'_'+'surveyconfig.cfg'

        survey_config_string = [toplevelstation,
                                'sampling = ' '{}'.format(sampling),
                                'acquired_by = ' '{}'.format(acquired_by),
                                'b_instrument_manufacturer = ' '{}'.format(b_instrument_manufacturer),
                                'b_instrument_model = ' '{}'.format(b_instrument_model),
                                'b_logger_type = ' '{}'.format(b_logger_type),
                                'b_instrument_type = ' '{}'.format(b_instrument_type),
                                'b_logger_gain = ' '{}'.format(b_logger_gain),
                                'b_instrument_amplification = ' '{}'.format(b_instrument_amplification),
                                'b_xaxis_azimuth = ' '{}'.format(b_xaxis_azimuth),
                                'b_yaxis_azimuth = ' '{}'.format(b_yaxis_azimuth),
                                'e_logger_type = ' '{}'.format(e_logger_type),
                                'e_logger_gain = ' '{}'.format(e_logger_gain),
                                'e_instrument_type = ' '{}'.format(e_instrument_type),
                                'e_instrument_amplification = ' '{}'.format(e_instrument_amplification),                        
                                'e_xaxis_azimuth = ' '{}'.format(e_xaxis_azimuth),
                                'e_xaxis_length = ' '{}'.format(e_xaxis_length),
                                'e_yaxis_azimuth = ' '{}'.format(e_yaxis_azimuth),
                                'e_yaxis_length = ' '{}'.format(e_yaxis_length),                                                
                                'data_logger_manufacturer = ' '{}'.format(data_logger_manufacturer),
                                'data_logger_type = ' '{}'.format(data_logger_type),
                                'lat = ' '{}'.format(lat),
                                'lon = ' '{}'.format(lon),
                                'elevation = ' '{}'.format(elevation),
                                'location = ' '{}'.format(location),
                                'station = ' '{}'.format(station),
                                'station_type = ' '{}'.format(station_type),
                                'date = ' '{}'.format(date),
                                'network = ' '{}'.format(network),
                                'hx = ' '{}'.format(hx),
                                'hy = ' '{}'.format(hy),
                                'hz = ' '{}'.format(hz),
                                'rr_box = ' '{}'.format(rr_box),
                                'rr_date = ' '{}'.format(rr_date),
                                'rr_hx = ' '{}'.format(rr_hx),
                                'rr_hy = ' '{}'.format(rr_hy),
                                'rr_lat = ' '{}'.format(rr_lat),
                                'rr_lon = ' '{}'.format(rr_lon),
                                'rr_station = ' '{}'.format(rr_station)]


        with open(survey_config_file_name,'w') as fp:
            for item in survey_config_string:
                fp.write("%s\n" % item)
            print('The survey config file has been written')
        
        
        ### run BIRRP with BIRRP script
        
        os.system('{} < {}'.format(birrp_location, birrp_script_name))

        print('birrp processing complete for: ' '{}_RR_{}'.format(station1,station2))
        
        ### remove local MTH5 files
        
        if os.path.exists(file1):
            os.remove(file1)
            print('removed' '{}'.format(file1))

        if os.path.exists(file2):
            os.remove(file2)
            print('removed ' '{}'.format(file2))
        
        ### create EDI file from BIRRP output j file, survey_config_file and birrp_config_file
              
        j2edi_obj = birrp.J2Edi()
        j2edi_obj.write_edi_file(station=station1, birrp_dir = new_directory, survey_config_fn = survey_config_file_name, birrp_config_fn = birrp_config_file_name, copy_path=None)
        edi_file = new_directory + '/' + station1 + '.edi'
        shutil.copy(edi_file,final_edi_dir)


        ### clean up processing directory
        ### NOTE: comment the following lines of code if you would 
        ### like to keep the intermediate BIRRP files for further analysis
        edi_files = glob.glob('*.edi')
        ffts = glob.glob('fft.*')
        coherences = glob.glob('*.2c2')
        rfs = glob.glob('*.rf')
        rps = glob.glob('*.rp')
        diags = glob.glob('*.diag')
        covariances = glob.glob('*.cov')
        jfiles = glob.glob('*.j')

        for edi_file in edi_files:
            os.remove(edi_file)
        for fft in ffts:
            os.remove(fft)        
        for coh in coherences:
            os.remove(coh)
        for rf in rfs:
            os.remove(rf)
        for rp in rps:
            os.remove(rp)
        for diag in diags:
            os.remove(diag)
        for cov in covariances:
            os.remove(cov)
        for jfile in jfiles:
            os.remove(jfile)


def rm_copyright_from_EDI_file(edi_dir):
    #### small script to remove the auto generated copyright message that MTpy produces when creating an EDI file
    os.chdir(edi_dir)
    line_number_to_delete = 58
    EDI_files = sorted(glob.glob('*.edi'))
    for EDI_file in EDI_files:
        print(EDI_file)
        with open(EDI_file,'r') as ef:
            contents = list(ef)
            #print(contents)
        del contents[line_number_to_delete - 1]
        with open(EDI_file,'w') as ef:
             for edi_line in contents:
                ef.write(edi_line)


### run birrp_processing on all station_combinations
birrp_processing(station_combinations,bp_workdir,bp_outdir)

comm.Barrier()


if rank==0:
    ### run rm_copyright_from_EDI_file on all produced EDI files
    rm_copyright_from_EDI_file(final_edi_dir)
    ### print total time to run script
    print('The script took {0} seconds !'.format(time.time()-startTime))


