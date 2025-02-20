#!/bin/bash

# Directory containing your raster files
input_dir="/Users/jr555/Documents/research/hydrology/drone/20160705/Raw"
output_dir="/Users/jr555/Documents/research/hydrology/drone/20160705/Reprojected"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .tif file in the input directory
for file in "$input_dir"/*.tif; do
    # Extract the base name of the file (without extension)
    base_name=$(basename "$file" .tif)
    
    # Run gdalwarp to reproject the file
    gdalwarp -t_srs EPSG:3413 "$file" "$output_dir/${base_name}_reprojected.tif"
    
    echo "Processed $file"
done