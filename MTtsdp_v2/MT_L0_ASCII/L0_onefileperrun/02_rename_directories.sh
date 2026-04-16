#!/bin/bash

TARGET_DIR="/g/data/<xy12>/<abc123>/AusLAMP/South_Australia/SE_SA_survey/level_0/Concatenated_Time_Series_ASCII_per_day"

# Loop through directories starting with "SA0"
for dir in "$TARGET_DIR"/SA0*; do
  # Check if it is a directory
  if [ -d "$dir" ]; then
    # Extract the base name
    base=$(basename "$dir")
    # Rename by removing the "0" after "SA"
    new_name="${base/SA0/SA}"
    # Perform the rename
    mv "$dir" "$TARGET_DIR/$new_name"
    echo "Renamed '$base' to '$new_name'"
  fi
done
