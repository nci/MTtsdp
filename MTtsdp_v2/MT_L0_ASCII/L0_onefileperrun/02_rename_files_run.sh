#!/bin/bash

# Define the parent directory containing all target files
PARENT_DIR="/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/Maralinga_Far_West_Coast_survey_LEMI/level_0/Concatenated_Time_Series_ASCII_per_run/SA53-2A"

# Set the station code prefix (e.g., 'SA63-2')
STATION_PREFIX="SA53_2A"

# Set the original filename pattern to look for (e.g., 'SA063-2')
ORIGINAL_PATTERN="SA53-2A"

# --- End of Configuration ---

# Find all files within PARENT_DIR (including subdirectories) that start with ORIGINAL_PATTERN
find "$PARENT_DIR" -type f -name "${ORIGINAL_PATTERN}*" | while read -r file; do
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
