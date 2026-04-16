#!/usr/bin/env python3
import sys
import os

MT_L0_ASCII = os.environ.get("MT_L0_ASCII")
supplements = os.environ.get("supplements")
survey_name = os.environ.get("survey_name")

# print(type(MT_L0_ASCII))
# print(type(supplements))
# print(type(survey_name))

# print(MT_L0_ASCII)
# print(supplements)
# print(survey_name)

log_file = supplements + "/" + survey_name + "/log_human.log"  # Define human-readable log file
log_file_machine = supplements + "/" + survey_name + "/log_machine.log"  # Define machine-readable log file

# Ensure log files are replaced instead of appended
open(log_file, "w").close()
open(log_file_machine, "w").close()

def log_print(out, pre=0, msg='empty message'):
    if out:
        with open(log_file, "a") as f:
            f.write('\t'*pre + msg + "\n")

def log_print2(size, path):
    path_parts = path.split('/')
    if len(path_parts) > 2:
        station = path_parts[1]  # Extract the station (first element after root)
        day_number = path_parts[2]  # Extract the number after the second '/'
    else:
        station = "UNKNOWN"
        day_number = "UNKNOWN"
    with open(log_file_machine, "a") as f:
        f.write(f"{size},{station},{day_number}\n")

class check_wc:
    def __init__(self, fnm, debug=False):
        self.fnm = fnm
        self.debug = debug
        self.vol = dict()
        self.rd()
        self.wrong_sta = set()
    
    def rd(self):
        log_print(True, 0, '# Input data from %s ...' % self.fnm)
        for line in open(self.fnm):
            n_sample, data_fnm = line.strip().split()
            n_sample = int(n_sample)
            sta, day, tmp = data_fnm.split('/')[-3:]
            dat, com = tmp.split('.')
            if sta not in self.vol:
                self.vol[sta] = dict()
            if com not in self.vol[sta]:
                self.vol[sta][com] = dict()
            self.vol[sta][com][day] = {'size': n_sample, 'nm': dat}
            self.rootdir = '/'+'/'.join(open(self.fnm).readline().split()[1].split('/')[:-3]) + '/'
    
    def get_path(self, com, sta, day):
        try:
            path = '/'.join([self.rootdir, sta, day, self.vol[sta][com][day]['nm']]) + '.' + com
            return path
        except:
            log_print(True, 0, 'Err: cannot find path for file %s %s %s' % (com, sta, day))
    
    def check_compare_com_at_same_station(self):
        log_print(True, 0, '# Checking matchness between different components for number of days...')
        for sta, v1 in self.vol.items():
            com_lst = sorted(list(v1.keys()))
            ndays_lst = [len(v1[com]) for com in com_lst]
            if ndays_lst.count(ndays_lst[0]) != len(ndays_lst):
                log_print(True, 0, 'Err: unmatched number of days between components. %s \n\t %s %s' % (sta, ', '.join(com_lst), ', '.join(['%d' % it for it in ndays_lst])))
                self.wrong_sta.add(sta)
    
    def check_wrong_day_name(self):
        log_print(True, 0, '# Checking wrong day names...')
        for sta, v1 in self.vol.items():
            for com, v2 in v1.items():
                for a_day in v2.keys():
                    if any(c not in '0123456789' for c in a_day):
                        path = self.get_path(com, sta, a_day).replace(self.rootdir, '')
                        log_print(True, 1, 'Err: wrong day name. %s' % path)
                        self.wrong_sta.add(sta)
                        break
    
    def check_gap_day(self):
        log_print(True, 0, '# Checking gap of days...')
        for sta in self.vol:
            for com in self.vol[sta]:
                day_lst = sorted(self.vol[sta][com].keys())
                try:
                    start = int(day_lst[0]) 
                    for idx, a_day in enumerate(day_lst):
                        if start + idx != int(a_day):
                            log_print(True, 1, 'Err: gap of days. %s %s; %s' % (sta, com, ', '.join(day_lst)))
                            self.wrong_sta.add(sta)
                except:
                    log_print(True, 1, 'Err: gap of days. %s %s; %s' % (sta, com, ', '.join(day_lst)))
                    self.wrong_sta.add(sta)
    
    def check_match_samples(self):
        log_print(True, 0, '# Checking number of samples per day')
        for sta in self.vol:
            for com in self.vol[sta]:
                day_lst = sorted(self.vol[sta][com].keys())
                for day in day_lst[1:-1]:
                    if self.vol[sta][com][day]['size'] not in [86400, 864000]:
                        path = self.get_path(com, sta, day).replace(self.rootdir, '')
                        log_print(True, 1, 'Warning: please check the number of samples %d %s' % (self.vol[sta][com][day]['size'], path))
                        log_print2(self.vol[sta][com][day]['size'], path)
                        self.wrong_sta.add(sta)
    
    def run(self):
        self.check_compare_com_at_same_station()
        self.check_wrong_day_name()
        self.check_gap_day()
        self.check_match_samples()
        log_print(True, 0, '\n\n#PLEASE check stations: %s' % (', '.join(self.wrong_sta)))

if __name__ == "__main__":
    app = check_wc(supplements + "/" + survey_name + "/list_fileperrun.log", debug=True)
    app.run()

