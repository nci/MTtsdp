#!/bin/bash

# Define the parent directory containing all target directories
PARENT_DIR="/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/SE_SA_survey/level_0/Concatenated_Time_Series_ASCII_per_day/SA03_2"


# Set the station code prefix (e.g., 'SA63-2')
STATION_PREFIX="SA03_2_"

# Set the original filename pattern to look for (e.g., 'SA063-2')
ORIGINAL_PATTERN="SA03_2"


# --- End of Configuration ---

# Loop through all subdirectories in the parent directory
find "$PARENT_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
  # Find files starting with the original pattern
  find "$dir" -type f -name "${ORIGINAL_PATTERN}*" | while read -r file; do
    # Extract directory and filename
    dir_path=$(dirname "$file")
    filename=$(basename "$file")
    
    # Replace the original pattern at the start of the filename with the station prefix
    new_filename=$(echo "$filename" | sed "s/^${ORIGINAL_PATTERN}/${STATION_PREFIX}/")
    
    # Rename the file if the filename has changed
    if [ "$filename" != "$new_filename" ]; then
      mv "$file" "$dir_path/$new_filename"
      echo "Renamed: $filename -> $new_filename"
    fi
  done
done
