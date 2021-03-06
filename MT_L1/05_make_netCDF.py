#!/usr/bin/env python3

import sys
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import pytz, datetime
import calendar
from dateutil.parser import parse

def mpi_log_print(out, pre=0, msg='empty message', file= sys.stdout, flush=True):
    if out:
        print('\t'*pre + msg, file= file, flush=flush)
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
class parameter_set:
    def __init__(self, key_fnm = '00_raw_ASCII_new/WA_and_SA_metadata/WA/auslamp_metadata_template.txt',
                    sta_fnm_template = '00_raw_ASCII_new/WA_and_SA_metadata/WA/%s.txt'):
        self.key_fnm = key_fnm
        self.sta_template= sta_fnm_template
        self.keys = [it.split('=')[0] for it in open(self.key_fnm, 'r')]
        self.vol =dict()
    def run(self, sta_set):
        self.sta_set = sta_set
        for sta in self.sta_set:
            if sta not in self.vol:
                self.vol[sta] = dict()
            fnm = self.sta_template % (sta)
            for k, v in zip(self.keys, open(fnm, 'r') ):
                k = k.strip()
                #print(k, v)
                self.vol[sta][k] = v.strip()
    def get_value(self, sta, key, raw=True):
        """
        Return the raw string of the value
        """
        try:
            v = self.vol[sta][key]
            if not raw:
                return float(v)
            else:
                return v
        except:
            mpi_log_print(True, 0, '`NA` Warning: cannot get values. %s %s ' % (sta, key), file=sys.stderr, flush=True )
            raise Exception
            
class make_nc:
    def __init__(self, timezone, fnm, para, channel_template, tempature_template, outname_template, lognm='', exclude_sta = set() ):
        """
        fnm: Filename of the *.log generated by `wc -l *` that is described in README.
        exclude_sta: stations to be excluded.
        """
        self.comm = MPI.COMM_WORLD
        self.mpi_size = self.comm.Get_size()
        self.rank     = self.comm.Get_rank()
        #self.rank = 0
        self.mpi_log = open('%s_%03d.log' % (lognm, self.rank), 'w' )
        #
        self.timezone = timezone
        #
        self.fnm = fnm
        self.exclude_sta = exclude_sta
        #
        self.para = para
        #
        self.outname_template = outname_template
        #
        self.channel_template = channel_template
        self.tempature_template = tempature_template
        #
        self.sta_vol = dict()
        #
        self.get_sta_set()
        #
    def run(self, cut_head=1, cut_end=1):
        self.para.run(self.sta_vol)
        stalst = sorted( self.sta_vol.keys() )
        chunk_size = len(stalst) // self.mpi_size + 1
        idx1 = chunk_size* self.rank
        idx2 = idx1 + chunk_size
        mpi_log_print(True, 0, '# Run nc ... %s' % (' '.join(stalst[idx1: idx2]) ) , file=self.mpi_log, flush=True )
        for sta in stalst[idx1: idx2]:
            self.run_single_sta(sta, cut_head, cut_end)
    def get_sta_set(self):
        """
        Get a set of station from `self.fnm`.
        """
        mpi_log_print(True, 0, '# Input data from %s ...' % self.fnm, file=self.mpi_log, flush=True )
        mpi_log_print(True, 1, 'ignore stations: %s' % (', '.join(list(self.exclude_sta))) , file=self.mpi_log, flush=True )
        for line in open(self.fnm):
            n_sample, data_fnm = line.strip().split()
            n_sample = int(n_sample)
            sta, day, tmp = data_fnm.split('/')[-3:]
            timestamp = tmp.split('_')[-1].split('.')[0]
            if sta in self.exclude_sta:
                continue # get rid of excluded stations
            if sta not in self.sta_vol:
                self.sta_vol[sta] = dict()
            self.sta_vol[sta][day] = timestamp
    def run_single_sta(self, stanm, cut_head=1, cut_end=1):
        """
        """
        mpi_log_print(True, 1, 'Run nc for sta: %s ...' % stanm, file=self.mpi_log, flush=True )
        daylst = sorted(self.sta_vol[stanm].keys() )[cut_head:-cut_head]
        mpi_log_print(True, 2, 'days(%d): %s ==> %s' % (len(daylst), daylst[0], daylst[-1] ), file=self.mpi_log, flush=True )
        out_name = self.outname_template % (stanm, int(daylst[0]), int(daylst[-1])  )
        dataset_P = Dataset(out_name, 'w', format='NETCDF4')
        ### global info
        dataset_P.info = 'All parameters are of raw values from metadata, whereas the time series data is processed by cut (%d %d), merge, rotation.' % (cut_head, cut_end)
        self.make_variable_global(stanm, dataset_P)### ch EX
        ### ch
        for ch, func in zip(
            ['ex', 'ey', 'bx', 'by', 'bz'],
            [self.make_variable_ex, self.make_variable_ey, self.make_variable_bx, self.make_variable_by, self.make_variable_bz]
        ):
            bin_fnm = self.channel_template % (stanm, cut_head, cut_end, ch.upper() )
            mpi_log_print(True, 2, bin_fnm, file=self.mpi_log, flush=True )
            arr = np.fromfile(bin_fnm, dtype='float32')
            d = dataset_P.createDimension(ch, arr.size)
            dd = dataset_P.createVariable(ch, np.float32, (ch,) )
            func(stanm, dd)
            dd[:] = arr
        ### time
        unixtime = timestamp_to_timeaxis(self.sta_vol[stanm][daylst[0] ], self.timezone)
        t1 = float(unixtime)
        number_of_samples = arr.size
        seconds = number_of_samples / 1.0 # 1Hz here after downsampling
        t2 = t1 + seconds
        time1 = np.linspace(t1, t2, num = number_of_samples)
        #
        time = dataset_P.createDimension('time', len(time1))
        time = dataset_P.createVariable('time',np.float64,('time',))
        time.units = 'seconds'
        time.long_name = 'time'
        time.standard_name = 'time'
        try:
            time.sampling_rate_Hz = '1' 
        except Exception:
            pass
        time[:] = time1
    def make_variable_ex(self, stanm, ex):
        ex.units = 'mV/m'
        ex.long_name = 'electric field in N-S orientation'
        ex.standard_name = 'electric field'
        ##
        try:
            ex.sampling_rate_Hz                         = '1'
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_NS_electrode_dipole_length_m')
            if  tmp_var != 'NA':
                ex.NS_electrode_dipole_length_m             = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_NS_electrode_azimuth_dec_deg')
            if  tmp_var != 'NA':
                ex.NS_electrode_azimuth_dec_deg             = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_NS_field_gain')
            if tmp_var != 'NA':
                ex.NS_field_gain = tmp_var                        
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_North_length_m')
            if tmp_var != 'NA':
                ex.North_length_m                           = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_North_electrode_number')
            if tmp_var != 'NA':
                ex.North_electrode_number                   = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_North_resistance_to_ground_kOhms')
            if tmp_var != 'NA':
                ex.North_resistance_to_ground_kOhms         = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_South_length_m')
            if tmp_var != 'NA':
                ex.South_length_m                           = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_South_electrode_number')
            if tmp_var != 'NA':
                ex.South_electrode_number                   = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_South_resistance_to_ground_kOhms')
            if tmp_var != 'NA':
                ex.South_resistance_to_ground_kOhms         = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'ground_center_electrode_number')
            if tmp_var != 'NA':
                ex.ground_center_electrode_number           = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_N_S_center_resistance_to_ground_kOhms')
            if tmp_var != 'NA':
                ex.N_S_center_resistance_to_ground_kOhms    = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_NS_AC_mV')
            if tmp_var != 'NA':
                ex.NS_AC_mV                                 = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'Ex_NS_DC_mV')
            if tmp_var != 'NA':
                ex.NS_DC_mV                                 = tmp_var
        except Exception:
            pass
        try:    
            tmp_var = self.para.get_value(stanm, 'MT_recorder_Ex_NS_voltage_mV')
            if tmp_var != 'NA':
                ex.MT_recorder_Ex_NS_voltage_mV             = tmp_var
        except Exception:
            pass
    def make_variable_ey(self, stanm, ey):
        ey.units = 'mV/m'
        ey.long_name = 'electric field in E-W orientation'
        ey.standard_name = 'electric field'
        #
        try:
            ey.sampling_rate_Hz                         = '1'
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_EW_electrode_dipole_length_m' )
            if tmp_var != 'NA':
                ey.EW_electrode_dipole_length_m             = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_EW_electrode_azimuth_dec_deg' )
            if tmp_var != 'NA':
                ey.EW_electrode_azimuth_dec_deg             = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_EW_field_gain' )
            if tmp_var != 'NA':
                ey.EW_field_gain                            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_East_length_m' )
            if tmp_var != 'NA':
                ey.East_length_m                            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_East_electrode_number' )
            if tmp_var != 'NA':
                ey.East_electrode_number                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var =  self.para.get_value(stanm, 'Ey_East_resistance_to_ground_kOhms' )
            if tmp_var != 'NA':
                ey.East_resistance_to_ground_kOhms          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_West_length_m' )
            if tmp_var != 'NA':
                ey.West_length_m                            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_West_electrode_number' )
            if tmp_var != 'NA':
                ey.West_electrode_number                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_West_resistance_to_ground_kOhms' )
            if tmp_var != 'NA':
                ey.West_resistance_to_ground_kOhms          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'ground_center_electrode_number' )
            if tmp_var != 'NA':
                ey.ground_center_electrode_number           = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_E_W_center_resistance_to_ground_kOhms' )
            if tmp_var != 'NA':
                ey.E_W_center_resistance_to_ground_kOhms    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_EW_AC_mV' )
            if tmp_var != 'NA':
                ey.EW_AC_mV                                 = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Ey_EW_DC_mV' )
            if tmp_var != 'NA':
                ey.EW_DC_mV                                 = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'MT_recorder_Ey_EW_voltage_mV' )
            if tmp_var != 'NA':
                ey.MT_recorder_Ey_EW_voltage_mV             = tmp_var
        except Exception:
            pass
    def make_variable_bx(self, stanm, bx):
        bx.units = 'nT'
        bx.long_name = 'magnetic field in X direction'
        bx.standard_name = 'magnetic field'
        try:
            bx.sampling_rate_Hz = '1'
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Bx_azimuth_degrees')
            if tmp_var != 'NA':
                bx.azimuth_degrees = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Bx_coil_number')
            if tmp_var != 'NA':
                bx.coil_number = tmp_var
        except Exception:
            pass
    def make_variable_by(self, stanm, by):
        by.units = 'nT'
        by.long_name = 'magnetic field in Y direction'
        by.standard_name = 'magnetic field'
        try:
            by.sampling_rate_Hz = '1'
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'By_azimuth_degrees')
            if tmp_var != 'NA':
                by.azimuth_degrees = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'By_coil_number')
            if tmp_var != 'NA':
                by.coil_number = tmp_var
        except Exception:
            pass
    def make_variable_bz(self, stanm, bz):
        bz.units = 'nT'
        bz.long_name = 'magnetic field in Z direction'
        bz.standard_name = 'magnetic field'
        try:
            bz.sampling_rate_Hz = '1'
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Bz_azimuth_degrees')
            if tmp_var != 'NA':
                bz.azimuth_degrees = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'Bz_coil_number')
            if tmp_var != 'NA':
                bz.coil_number = tmp_var
        except Exception:
            pass
    def make_variable_global(self, stanm, dataset_P):
        try:
            tmp_var = self.para.get_value(stanm, 'deployment_survey_prospect_name')
            if tmp_var != 'NA':
                dataset_P.title         = tmp_var
        except Exception:
            pass
        try:
            dataset_P.summary       = 'Long period time series for site %s of the %s survey' % (
                    self.para.get_value(stanm, 'site_number'), self.para.get_value(stanm, 'deployment_survey_prospect_name') )
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'operator')
            if tmp_var != 'NA':
                dataset_P.source        = tmp_var
        except Exception:
            pass
        try:
            dataset_P.date_created  = '@'+ self.timezone+ '@' + datetime.datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
        except Exception:
            pass
        try:
            dataset_P.Conventions   = 'ACDD-1.3'
        except Exception:
            pass                           
        try:
            dataset_P.license       = 'Creative Commons Attribution 4.0 International (CC BY 4.0)'
        except Exception:
            pass
        ### Deployment metadata
        try:
            tmp_var = self.para.get_value(stanm, 'site_number')
            if tmp_var != 'NA':
                dataset_P.site_number                                   = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'site_name')
            if tmp_var != 'NA':
                dataset_P.site_name                                     = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'location')
            if tmp_var != 'NA':
                dataset_P.location                                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'deployment_survey_prospect_name')
            if tmp_var != 'NA':
                dataset_P.deployment_survey_prospect_name               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'institution_agency')
            if tmp_var != 'NA':
                dataset_P.institution_agency                            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'AUSGRID_number')
            if tmp_var != 'NA':
                dataset_P.AUSGRID_number                                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'operator')
            if tmp_var != 'NA':
                dataset_P.operator                                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'recorder_latitude_dec_deg')
            if tmp_var != 'NA':
                dataset_P.recorder_latitude_dec_deg                     = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'recorder_longitude_dec_deg')
            if tmp_var != 'NA':
                dataset_P.recorder_longitude_dec_deg                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'recorder_GPS_projection_system')
            if tmp_var != 'NA':
                dataset_P.recorder_GPS_projection_system                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'field_GPS_latitude_dec_deg')
            if tmp_var != 'NA':
                dataset_P.field_GPS_latitude_dec_deg                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'field_GPS_longitude_dec_deg')
            if tmp_var != 'NA':
                dataset_P.field_GPS_longitude_dec_deg                   = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'GPS_easting_AMG')
            if tmp_var != 'NA':
                dataset_P.GPS_easting_AMG                               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'GPS_northing_AMG')
            if tmp_var != 'NA':
                dataset_P.GPS_northing_AMG                              = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'AMG_zone')
            if tmp_var != 'NA':
                dataset_P.AMP_zone                                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'GPS_projection_system')
            if tmp_var != 'NA':
                dataset_P.GPS_projection_system                         = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'altitude_or_elevation_from_GPS_or_instrument')
            if tmp_var != 'NA':
                dataset_P.altitude_or_elevation_from_GPS_or_instrument  = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'deployment_date')
            if tmp_var != 'NA':
                dataset_P.deployment_date                               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'LOCAL_start_hour')
            if tmp_var != 'NA':
                dataset_P.LOCAL_start_hour                              = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'LOCAL_start_minute')
            if tmp_var != 'NA':
                dataset_P.LOCAL_start_minute                            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'LOCAL_start_second')
            if tmp_var != 'NA':
                dataset_P.LOCAL_start_second                            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_date_from_instrument')
            if tmp_var != 'NA':
                dataset_P.UTC_date_from_instrument                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_year')
            if tmp_var != 'NA':
                dataset_P.UTC_start_year                                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_month')
            if tmp_var != 'NA':
                dataset_P.UTC_start_month                               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_day')
            if tmp_var != 'NA':
                dataset_P.UTC_start_day                                 = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_time')
            if tmp_var != 'NA':
                dataset_P.UTC_start_time                                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_time_hour')
            if tmp_var != 'NA':
                dataset_P.UTC_start_time_hour                           = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_time_minute')
            if tmp_var != 'NA':
                dataset_P.UTC_start_time_minute                         = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'UTC_start_time_second')
            if tmp_var != 'NA':
                dataset_P.UTC_start_time_second                         = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'deployment_Julian_day')
            if tmp_var != 'NA':
                dataset_P.deployment_Julian_day                         = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'recording_method')
            if tmp_var != 'NA':
                dataset_P.recording_method                              = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'MT_recorder_type_model')
            if tmp_var != 'NA':
                dataset_P.MT_recorder_type_model                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'magnetometer_type_model')
            if tmp_var != 'NA':
                dataset_P.magnetometer_type_model                       = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'electrode_type_model')
            if tmp_var != 'NA':
                dataset_P.electrode_type_model                          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'power_source_type_or_model')
            if tmp_var != 'NA':
                dataset_P.power_source_type_or_model                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'data_confidentiality')
            if tmp_var != 'NA':
                dataset_P.data_confidentiality                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'north_reference')
            if tmp_var != 'NA':
                dataset_P.north_reference                           = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'drift_calculation_end_time_GPS')
            if tmp_var != 'NA':
                dataset_P.drift_calculation_end_time_GPS            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'drift_calculation_end_time_instrument')
            if tmp_var != 'NA':
                dataset_P.drift_calculation_end_time_instrument     = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'instrument_drift_in_seconds')
            if tmp_var != 'NA':
                dataset_P.instrument_drift_in_seconds               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'battery_voltage_no_load_V')
            if tmp_var != 'NA':
                dataset_P.battery_voltage_no_load_V                 = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'battery_voltage_on_load_V')
            if tmp_var != 'NA':
                dataset_P.battery_voltage_on_load_V                 = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'battery_voltage_MT_data_logger_V')
            if tmp_var != 'NA':
                dataset_P.battery_voltage_MT_data_logger_V          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'MT_box_case_number')
            if tmp_var != 'NA':
                dataset_P.MT_box_case_number                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'MT_recorder_data_logger_number')
            if tmp_var != 'NA':
                dataset_P.MT_recorder_data_logger_number            = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'MT_recorder_interface_box_number')
            if tmp_var != 'NA':
                dataset_P.MT_recorder_interface_box_number          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'data_storage_device_number')
            if tmp_var != 'NA':
                dataset_P.data_storage_device_number                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'magnetometer_type')
            if tmp_var != 'NA':
                dataset_P.magnetometer_type                         = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'data_file_length_minutes_per_file')
            if tmp_var != 'NA':
                dataset_P.data_file_length_minutes_per_file         = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'NS_EW_resistance_kOhms')
            if tmp_var != 'NA':
                dataset_P.NS_EW_resistance_kOhms                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'NS_EW_voltage_mV')
            if tmp_var != 'NA':
                dataset_P.NS_EW_voltage_mV                          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'deployment_comments')
            if tmp_var != 'NA':
                dataset_P.deployment_comments                       = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'photo_taken_deployment')
            if tmp_var != 'NA':
                dataset_P.photo_taken_deployment                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'photo_identification_numbers_deployment')
            if tmp_var != 'NA':
                dataset_P.photo_identification_numbers_deployment   = tmp_var
        except Exception:
            pass
        #dataset_P.deployment_sheet_scan = deployment_sheet_scan
        try:
            tmp_var = self.para.get_value(stanm, 'calculated_date_for_retrieval')
            if tmp_var != 'NA':
                dataset_P.calculated_date_for_retrieval             = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_date_local')
            if tmp_var != 'NA':
                dataset_P.retrieval_date_local                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_local_start_year')
            if tmp_var != 'NA':
                dataset_P.retrieval_local_start_year                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_local_start_month')
            if tmp_var != 'NA':
                dataset_P.retrieval_local_start_month               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_local_start_day')
            if tmp_var != 'NA':
                dataset_P.retrieval_local_start_day                 = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_time_local')
            if tmp_var != 'NA':
                dataset_P.retrieval_time_local                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_local_start_hour')
            if tmp_var != 'NA':
                dataset_P.retrieval_local_start_hour                = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_local_start_minute')
            if tmp_var != 'NA':
                dataset_P.retrieval_local_start_minute              = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_local_start_second')
            if tmp_var != 'NA':
                dataset_P.retrieval_local_start_second              = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_date_UTC')
            if tmp_var != 'NA':
                dataset_P.retrieval_date_UTC                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_time_UTC')
            if tmp_var != 'NA':
                dataset_P.retrieval_time_UTC                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_Julian_day')
            if tmp_var != 'NA':
                dataset_P.retrieval_Julian_day                      = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_Ex_N_S_Center_resistance_kOhms')
            if tmp_var != 'NA':
                dataset_P.retrieval_Ex_N_S_Center_resistance_kOhms  = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_Ex_AC_mV')
            if tmp_var != 'NA':
                dataset_P.retrieval_Ex_AC_mV                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_Ex_DC_mV')
            if tmp_var != 'NA':
                dataset_P.retrieval_Ex_DC_mV                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_Ey_E_W_Center_resistance_kOhms')
            if tmp_var != 'NA':
                dataset_P.retrieval_Ey_E_W_Center_resistance_kOhms  = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_ey_AC_mV')
            if tmp_var != 'NA':
                dataset_P.retrieval_ey_AC_mV                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_Ey_DC_mV')
            if tmp_var != 'NA':
                dataset_P.retrieval_Ey_DC_mV                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_NS_EW_resistance_kOhms')
            if tmp_var != 'NA':
                dataset_P.retrieval_NS_EW_resistance_kOhms          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'MT_unit_operational')
            if tmp_var != 'NA':
                dataset_P.MT_unit_operational                       = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_battery_voltage_V')
            if tmp_var != 'NA':
                dataset_P.retrieval_battery_voltage_V               = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'electrode_cable_status')
            if tmp_var != 'NA':
                dataset_P.electrode_cable_status                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'mag_cable_status')
            if tmp_var != 'NA':
                dataset_P.mag_cable_status                          = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'retrieval_comments')
            if tmp_var != 'NA':
                dataset_P.retrieval_comments                        = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'days_operated')
            if tmp_var != 'NA':
                dataset_P.days_operated                             = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'photos_taken_retrieval')
            if tmp_var != 'NA':
                dataset_P.photos_taken_retrieval                    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'photo_identification_numbers_retrieval')
            if tmp_var != 'NA':
                dataset_P.photo_identification_numbers_retrieval    = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'site_data_quality_retrieval')
            if tmp_var != 'NA':
                dataset_P.site_data_quality_retrieval   = tmp_var
        except Exception:
            pass
        try:
            tmp_var = self.para.get_value(stanm, 'redeployed')
            if tmp_var != 'NA':
                dataset_P.redeployed                    = tmp_var
        except Exception:
            pass

if __name__ == "__main__":
    WA_sta_exclude_set = set( (
        'WA11','WA12','WA28', 'WASA352', # different machine
        'WA26', 'WASA302' # warnings reported by 01_check.py
    ) )

    SA_sta_exclude_set = set( (
            'SA246', 'SA299', 'SA351', 'SA320-2', 'SA227', 'SA347', 'SA247', 
            'SA250', 'SA324-2', 'SA326S', 'SA277', 'SA344-2'
        ) )
    #########
    timezone = "Australia/Perth"
    ######### WA
    para = parameter_set(key_fnm = '00_raw_ASCII_new/WA_and_SA_metadata/WA/auslamp_metadata_template.txt',
                            sta_fnm_template = '00_raw_ASCII_new/WA_and_SA_metadata/WA/%s.txt')
    app = make_nc(timezone,
        '01_workspace/WA.log', para, 
        '04_workspace/rotated_data_bin_WA/merged-%s-cut_%03d_%03d.%s.bin', 
        '03_workspace/merged_data_bin_WA/merged-%s-cut_%03d_%03d.ambientTemperature.bin',
        '05_workspace/WA/%s_%03d_%03d.nc',
        '05_workspace/WA_log_', WA_sta_exclude_set)
    app.run(cut_head=1, cut_end=1)
    ######### SA
    para = parameter_set(key_fnm = '00_raw_ASCII_new/WA_and_SA_metadata/SA/auslamp_metadata_template.txt',
                            sta_fnm_template = '00_raw_ASCII_new/WA_and_SA_metadata/SA/%s.txt')
    app = make_nc(timezone,
        '01_workspace/SA.log', para, 
        '04_workspace/rotated_data_bin_SA/merged-%s-cut_%03d_%03d.%s.bin', 
        '03_workspace/merged_data_bin_SA/merged-%s-cut_%03d_%03d.ambientTemperature.bin',
        '05_workspace/SA/%s_%03d_%03d.nc',
        '05_workspace/SA_log_', SA_sta_exclude_set)
    app.run(cut_head=1, cut_end=1)
