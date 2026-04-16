### USAGE:
### python 02_extract_translated_azimuth.py <survey_name>
### For example:
### python 02_extract_translated_azimuth.py Maralinga_Far_West_Coast_survey

import numpy as np
import json
import csv
import os
from mpi4py import MPI
import sys

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def read_ascii_file(filename):
    return np.loadtxt(filename)

# Check if survey name is provided
if len(sys.argv) < 2:
    print("Usage: python your_script.py <survey_name>")
    sys.exit(1)

survey = sys.argv[1]

level = 'level_0'
run_product = 'Concatenated_Time_Series_ASCII_per_run'
workdir_run = f'/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/{survey}/{level}/{run_product}'
output_dir = f'/g/data/<xy12>/<abc123>/supplementing_files/translated_azimuth/{survey}'
os.makedirs(output_dir, exist_ok=True)

# List all stations (done only by rank 0, then broadcast)
if rank == 0:
    stations = sorted([d for d in os.listdir(workdir_run) if os.path.isdir(os.path.join(workdir_run, d))])
else:
    stations = None

stations = comm.bcast(stations, root=0)

# Divide stations among processes
stations_per_proc = np.array_split(stations, size)
my_stations = stations_per_proc[rank]

# Each process processes its subset of stations
local_results = []
local_summary = {}

for station in my_stations:
    folder_run = os.path.join(workdir_run, station)
    bx_files = sorted([f for f in os.listdir(folder_run) if f.endswith('.BX')])
    by_files = sorted([f for f in os.listdir(folder_run) if f.endswith('.BY')])
    for bx_file, by_file in zip(bx_files, by_files):
        bx_filename = os.path.join(folder_run, bx_file)
        by_filename = os.path.join(folder_run, by_file)
        BX1 = read_ascii_file(bx_filename)
        BY1 = read_ascii_file(by_filename)
        meanBx = np.mean(BX1[5000:])
        meanBy = np.mean(BY1[5000:])
        z = meanBx + 1j*meanBy
        RotAngleB = np.angle(z)
        print(f"Process {rank} - Station: {station}, Azimuth: {RotAngleB:.4f}")
        # Store results locally
        local_results.append({'station': station, 'translated_azimuth': RotAngleB})
        # Store in local summary
        local_summary[station] = RotAngleB

# Gather all results at rank 0
all_results = comm.gather(local_results, root=0)
all_summary = comm.gather(local_summary, root=0)

if rank == 0:
    # Flatten the list of results
    combined_results = [item for sublist in all_results for item in sublist]
    # Merge summaries
    combined_summary = {}
    for d in all_summary:
        combined_summary.update(d)

    # Write CSV file
    csv_filename = 'station_azimuths.csv'
    csv_file = os.path.join(output_dir, csv_filename)
    with open(csv_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['station', 'translated_azimuth'])
        writer.writeheader()
        for row in combined_results:
            azimuth_str = f"{row['translated_azimuth']:.4f}"
            writer.writerow({'station': row['station'], 'translated_azimuth': azimuth_str})
    print(f"CSV file saved as {csv_file}")

    # Write JSON summary
    json_filename = 'station_azimuths.json'
    json_file = os.path.join(output_dir, json_filename)
    formatted_summary = {station: f"{azimuth:.4f}" for station, azimuth in combined_summary.items()}
    with open(json_file, 'w') as jsonfile:
        json.dump(formatted_summary, jsonfile, indent=4)
    print(f"JSON summary saved as {json_file}")
