#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute mean depth of meltwater ponds from Zhang et al. (2019)

"""

# Import packages
import rasterio
import numpy as np

# Define path
path = "/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/"

# Define functions
def compute_raster_mean(raster_path):
    with rasterio.open(raster_path) as src:
        mean_values = []
        for block_index, window in src.block_windows(1):  # Read in chunks
            data = src.read(1, window=window, masked=True)  # Read with masking
            if data.count() > 0:  # Ignore empty chunks
                mean_values.append(data.mean())
        
        return np.mean(mean_values) if mean_values else None  # Compute overall mean

def count_ones_in_raster(raster_path):
    with rasterio.open(raster_path) as src:
        count_ones = 0
        for block_index, window in src.block_windows(1):  # Read in chunks
            data = src.read(1, window=window, masked=True)  # Read with masking
            if data.count() > 0:  # Ignore empty chunks
                count_ones += (data == 1).sum()  # Count the number of 1s in the chunk
        return count_ones
    
# Define functions
def compute_raster_std(raster_path):
    with rasterio.open(raster_path) as src:
        std_values = []
        for block_index, window in src.block_windows(1):  # Read in chunks
            data = src.read(1, window=window, masked=True)  # Read with masking
            if data.count() > 0:  # Ignore empty chunks
                std_values.append(data.std())
        return np.mean(std_values) if std_values else None  # Compute overall mean
    
#%%

raster_path = path + '/zhang/surface_water_depth_2019.tif'
mean_value = compute_raster_mean(raster_path)
std_value = compute_raster_std(raster_path)
print(f"Mean raster value: {mean_value}")
print(f"Std. dev. raster value: {std_value}")

#%%

crevasse_path = path + '/zhang/water-filled_crevasses_2019.tif'
sum_crevasse = count_ones_in_raster(crevasse_path)

water_path = path + '/zhang/surface_water_mask_2019.tif'
sum_water = count_ones_in_raster(water_path)
print('Area of all water in 2019 from Zhang et al. is %.2f' %(sum_water/10000))
print('Area of water-filled crevasses in 2019 from Zhang et al. is %.2f' %(sum_crevasse/10000))
#%%


