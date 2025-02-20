#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute mean depth of meltwater ponds from Zhang et al. (2019)

"""


import rasterio
import numpy as np

def compute_raster_mean(raster_path):
    with rasterio.open(raster_path) as src:
        mean_values = []
        for block_index, window in src.block_windows(1):  # Read in chunks
            data = src.read(1, window=window, masked=True)  # Read with masking
            if data.count() > 0:  # Ignore empty chunks
                mean_values.append(data.mean())
        
        return np.mean(mean_values) if mean_values else None  # Compute overall mean

raster_path = "/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/surface_water_depth_2019.tif"
mean_value = compute_raster_mean(raster_path)
print(f"Mean raster value: {mean_value}")
