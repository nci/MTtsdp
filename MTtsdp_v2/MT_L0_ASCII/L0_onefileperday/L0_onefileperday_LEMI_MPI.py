import mth5
from datetime import datetime
import re
import os
import xarray as xr
import time
from mpi4py import MPI

startTime = time.time()

# Initialize MPI
comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

### Define final directory where outputs go
final_dir = '/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/Maralinga_Far_West_Coast_survey_LEMI/level_0/Concatenated_Time_Series_ASCII_per_day'

### Define the working directory where the LEMI data lies for each station
working_directory = '/g/data/<xy12>/<abc123>/raw_data/LEMI/Maralinga_LEMI'

station_dirs = [
    os.path.join(working_directory, dirname)
    for dirname in os.listdir(working_directory)
    if os.path.isdir(os.path.join(working_directory, dirname))
]

def get_julian_day(date_str):
    date_obj = datetime.strptime(date_str, "%Y%m%d")
    return date_obj.timetuple().tm_yday

def convert_date_string(date_str):
    year = date_str[0:4]
    month = date_str[4:6]
    day = date_str[6:8]
    hour = date_str[8:10]
    minute = date_str[10:12]
    dt = datetime(int(year), int(month), int(day), int(hour), int(minute))
    new_format = dt.strftime('%y%m%d%H%M')
    return new_format + '00'

def save_values(file_path, values):
    with open(file_path, 'w') as file:
        file.writelines(f"{value}\n" for value in values)

def extract_date_from_filename(filename):
    pattern = r'/(\d{8})\d{4}\.TXT$'
    match = re.search(pattern, filename)
    if match:
        return match.group(1)
    else:
        return None

def process_lemi_files(file_list, station, top_dir):
    for lemi_example_file in file_list:
        filename = os.path.basename(lemi_example_file)
        timestamp_str = os.path.splitext(filename)[0]

        # Read the LEMI file
        read_lemi = mth5.read_file(lemi_example_file, file_type='lemi424')
        dataset = read_lemi.dataset

        # extract the electromagnetic field data
        bx = dataset['bx']
        by = dataset['by']
        bz = dataset['bz']
        e1 = dataset['e1']
        e2 = dataset['e2']
        temperature_e = dataset['temperature_e']
        temperature_h = dataset['temperature_h']

        # extract the raw time series values
        raw_bx = bx.values
        raw_by = by.values
        raw_bz = bz.values
        raw_e1 = e1.values
        raw_e2 = e2.values
        raw_temp_e = temperature_e.values
        raw_temp_h = temperature_h.values

        # extract the date, earth data logger (EDL) date string and Julian day from filename
        date_str = extract_date_from_filename(lemi_example_file)
        edl_date_string = convert_date_string(timestamp_str)
        julian_day = str(get_julian_day(date_str))

        # construct directory path
        directory_path = os.path.join(top_dir, station, julian_day)
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
            if rank == 0:
                print(f"Directory created: {directory_path}")
        # Save values
        bx_filename = f"{station}_{edl_date_string}.BX"
        by_filename = f"{station}_{edl_date_string}.BY"
        bz_filename = f"{station}_{edl_date_string}.BZ"
        ex_filename = f"{station}_{edl_date_string}.EX"
        ey_filename = f"{station}_{edl_date_string}.EY"
        temp_e_filename = f"{station}_{edl_date_string}.TPE"
        temp_h_filename = f"{station}_{edl_date_string}.TPH"

        file_path_bx = os.path.join(directory_path, bx_filename)
        file_path_by = os.path.join(directory_path, by_filename)
        file_path_bz = os.path.join(directory_path, bz_filename)
        file_path_ex = os.path.join(directory_path, ex_filename)
        file_path_ey = os.path.join(directory_path, ey_filename)
        file_path_tpe = os.path.join(directory_path, temp_e_filename)
        file_path_tph = os.path.join(directory_path, temp_h_filename)

        save_values(file_path_bx, raw_bx)
        save_values(file_path_by, raw_by)
        save_values(file_path_bz, raw_bz)
        save_values(file_path_ex, raw_e1)
        save_values(file_path_ey, raw_e2)
        save_values(file_path_tpe, raw_temp_e)
        save_values(file_path_tph, raw_temp_h)


# Distribute station_dirs among ranks
# Use simple round-robin assignment
for i, station_dir in enumerate(station_dirs):
    if i % size == rank:
        station = os.path.basename(station_dir)
        file_list = sorted([
            os.path.join(station_dir, filename)
            for filename in os.listdir(station_dir)
            if filename.endswith('.TXT') and os.path.isfile(os.path.join(station_dir, filename))
        ])
        # Proceed only if there are files
        if file_list:
            process_lemi_files(file_list, station, final_dir)

# Optional: synchronize all processes at the end
comm.Barrier()

print(f"Rank {rank} finished processing.")

# If you want to measure total runtime
if rank == 0:
    total_time = time.time() - startTime
    print(f"Total runtime: {total_time:.2f} seconds")
