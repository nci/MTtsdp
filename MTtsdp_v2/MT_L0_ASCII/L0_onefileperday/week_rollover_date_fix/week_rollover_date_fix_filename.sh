#!/bin/bash

# Define the list of folders
folders=("SA28-2B" "SA33-2A" "SA38-2B" "SA46-2B" "SA50-2A" "SA53-2B") 

for folder in "${folders[@]}"; do
    echo "Processing folder: $folder"

    # Change to the target directory
    if [ -d "$folder" ]; then
        cd "$folder" || { echo "Failed to enter directory '$folder'"; continue; }
    else
        echo "Directory '$folder' not found, skipping."
        continue
    fi

    # Loop over all .TXT files in the current folder
    for filename in *.TXT; do
        # Check if the file exists (in case no .TXT files are found)
        [ -e "$filename" ] || continue

        # Extract the date part (first 8 characters)
        date_part=${filename:0:8}

        # Convert to ISO format
        date_iso=$(date -d "${date_part}" +%Y-%m-%d) || { echo "Invalid date in filename '$filename'"; continue; }

        # Add 7168 days
        new_date=$(date -d "$date_iso +7168 days" +%Y%m%d)

        # Construct new filename
        new_filename="${new_date}0000.TXT"

        # Rename the file
        mv "$filename" "$new_filename"

        echo "Renamed '$filename' to '$new_filename'"
    done

    # Return to the original directory before processing the next folder
    cd - > /dev/null
done

