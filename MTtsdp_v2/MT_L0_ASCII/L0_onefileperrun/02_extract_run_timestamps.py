import sys
import os
from datetime import datetime, timezone, timedelta
import re
import json
import csv

save_csv = True


# Check if survey name is provided
if len(sys.argv) < 3:
    print("Usage: python your_script.py <survey_name> <instrument_type>")
    print("Where <instrument_type> is 'EDL' or 'LEMI'")
    sys.exit(1)

survey = sys.argv[1]
instrument_type = sys.argv[2]

# Validate instrument_type
if instrument_type not in ("EDL", "LEMI"):
    print("Error: instrument_type must be 'EDL' or 'LEMI'")
    sys.exit(1)


level = 'level_0'
run_product = 'Concatenated_Time_Series_ASCII_per_run'
day_product = 'Concatenated_Time_Series_ASCII_per_day'
workdir_run = f'/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/{survey}/{level}/{run_product}'
workdir_day = f'/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/{survey}/{level}/{day_product}'
output_dir = f'/g/data/<xy12>/<abc123>/supplementing_files/run_timestamps/{survey}'
os.makedirs(output_dir, exist_ok=True)

stations = sorted([d for d in os.listdir(workdir_run) if os.path.isdir(os.path.join(workdir_run, d))])

# Initialize a list to store all summaries
summary_list = []

for station in stations:
    print(f"\nProcessing station: {station}")
    folder_run = os.path.join(workdir_run, station)
    
    # Initialize lists to store start/end times per run
    all_start_times = []
    all_end_times = []
    
    # 1. Find all .BX files in the run directory
    bx_files = sorted([f for f in os.listdir(folder_run) if f.endswith('.BX')])
    
    # 2. For each .BX file, extract station, JD1, JD2
    run_info_list = []
    
    for filename in bx_files:
        match = re.match(r'^([A-Za-z0-9_]+)-(\d+)_(\d+)\.BX$', filename)
        if match:
            station = match.group(1)
            JD1 = int(match.group(2))
            JD2 = int(match.group(3))
            run_info_list.append({
                'filename': filename,
                'station': station,
                'JD1': JD1,
                'JD2': JD2
            })
        else:
            print(f'Filename {filename} does not match expected pattern.')
    
    # Process each run (JD1, JD2)
    for run_idx, run in enumerate(run_info_list,start=1):
        JD1 = run['JD1']
        JD2 = run['JD2']
        filename = run['filename']
        station = run['station']
        
        print("station=", station)
        print("filename=", filename)

        JD1_str = str(JD1).zfill(3)
        JD2_str = str(JD2).zfill(3)
        
        print("JD1=", JD1)
        print("JD2=", JD2)
        
        print("JD1_srt=", JD1_str)
        print("JD2_srt=", JD2_str)

        target_dir_JD1 = os.path.join(workdir_day, station, JD1_str)
        print("target_dir_JD1 =", target_dir_JD1) 

        #print(f"Navigating to directory: {target_dir_JD1}")
    
        try:
            files = os.listdir(target_dir_JD1)
        except FileNotFoundError:
            print(f"Directory not found: {target_dir_JD1}")
            sys.exit(1)    
    
        filename_JD1 = [f for f in files if f.endswith('.BX')]
        
        for file in filename_JD1:
            print(file)
            basename = os.path.basename(file)
            name_part, ext = os.path.splitext(basename)
            
            underscore_count = basename.count('_')

            print("underscore_count=",underscore_count)
            
            if underscore_count == 2:
                #first_underscore_index = basename.find('_')
                #second_underscore_index = basename.find('_', first_underscore_index + 1)
                #part1 = basename[:second_underscore_index]
                #part2 = basename[second_underscore_index + 1:]
                first_idx = name_part.find('_')
                second_idx = name_part.find('_', first_idx + 1)

                part1 = name_part[:second_idx]
                part2 = name_part[second_idx + 1:]
                
                parts = [part1, part2]
            else:
                parts = name_part.split('_')
                
            print("parts1=",parts)
            
            if len(parts) != 2:
                raise ValueError("Filename does not match expected format.")
            
            timestamp_str = parts[1]
            dt = datetime.strptime(timestamp_str, '%y%m%d%H%M%S')
            dt = dt.replace(tzinfo=timezone.utc)
            formatted_start_time = dt.isoformat(timespec='microseconds')
            all_start_times.append(formatted_start_time)
    
        target_dir_JD2 = os.path.join(workdir_day, station, JD2_str)  
        #print(f"Navigating to directory: {target_dir_JD2}")
    
        try:
            files_jd2 = os.listdir(target_dir_JD2)
        except FileNotFoundError:
            print(f"Directory not found: {target_dir_JD2}")
            sys.exit(1)
    
        filename_JD2 = [f for f in files_jd2 if f.endswith('.BX')]
        #print(f"filename_JD2:{filename_JD2}")
    
        for file in filename_JD2:
            filepath_JD2 = os.path.join(target_dir_JD2, file)
            #print(f"filepath_JD2:{filepath_JD2}")
            with open(filepath_JD2, 'r') as f:
                line_count = sum(1 for _ in f)
                if instrument_type == "EDL":
                    seconds_in_file = int(line_count / 10)
                elif instrument_type == "LEMI":
                    seconds_in_file = int(line_count)
            basename = os.path.basename(file)
            name_part, ext = os.path.splitext(basename)
            
            underscore_count = basename.count('_')
            
            print("underscore_count_2=",underscore_count)
            
            if underscore_count == 2:
                first_idx = name_part.find('_')
                second_idx = name_part.find('_', first_idx + 1)
                part1 = name_part[:second_idx]
                part2 = name_part[second_idx + 1:]
                parts = [part1, part2]
            else:
                parts = name_part.split('_')
            
            print("parts2=",parts)
            
            if len(parts) != 2:
                raise ValueError("Filename does not match expected format.")       
    
            timestamp_str = parts[1]
            dt = datetime.strptime(timestamp_str, '%y%m%d%H%M%S')
            dt = dt.replace(tzinfo=timezone.utc)
            #formatted_end_time_filename = dt.isoformat(timespec='microseconds')
            #end_time = datetime.fromisoformat(formatted_end_time_filename)
            end_time = dt
            
            if seconds_in_file == 86400:
                seconds_in_file -= 0.000001
            
            new_end_time = end_time + timedelta(seconds=seconds_in_file)
            formatted_end_time = new_end_time.isoformat()
            #print(new_end_time)
            all_end_times.append(formatted_end_time)

        # Store summary info for this run
        summary_list.append({
            'Station': station,
            'Run': run_idx,
            'end_time': formatted_end_time,
            'start_time': formatted_start_time
        })

# After processing all stations, write the summary to a JSON file
summary_filename = os.path.join(output_dir, 'stations_summary.json')
with open(summary_filename, 'w') as json_file:
    json.dump(summary_list, json_file, indent=4)

print(f"\nSummary JSON saved to: {summary_filename}")

# Save to CSV if flag is set
if save_csv:
    csv_filename = os.path.join(output_dir, 'stations_summary.csv')
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['Station', 'Run', 'End Time', 'Start Time'])
        # Write rows
        for entry in summary_list:
            writer.writerow([entry['Station'], entry['Run'], entry['end_time'], entry['start_time']])
    print(f"\nSummary CSV saved to: {csv_filename}")

# Optional: print the summary table
#print(f"\n{'Station':<10} {'Run':<5} {'end_time':<26} {'start_time':<26}")
#for entry in summary_list:
#    print(f"{entry['Station']:<10} {entry['Run']:<5} {entry['end_time']:<26} {entry['start_time']:<26}")
