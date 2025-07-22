#!/bin/bash

# Directory containing your input raster files
input_dir="/Users/jr555/Documents/research/hydrology/drone/russell/Reprojected"
output_dir="/Users/jr555/Documents/research/hydrology/drone/russell/Orthomosaics"

# Loop through each .tif file in the input directory
for file in "$input_dir"/*.tif; do
    # Run gdal_retile.py for each file
    gdal_retile.py -targetDir "$output_dir" -ps 2000 2000 -co "TILED=YES" "$file"
    
    echo "Processed $file"
done