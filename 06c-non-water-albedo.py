#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Quantify local albedo variability not due to surface water witin 6.25km2 area

"""


# Import packages
import rioxarray as rio
import xarray as xr
import numpy as np
import pandas as pd
from scipy.signal import convolve2d

#%%

# Define year
year = str(2019)

# Define sample size
sample_size = 100

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Read raster data
sw = rio.open_rasterio(path1 + 'zhang/surface_water_mask_' + year + '_500m.tif')

# Replace NaN values with zeros
sw = sw.fillna(0)

# Read MODIS product
mcd = rio.open_rasterio(path1 + 'mcd-' + year + '.tif')
mod = rio.open_rasterio(path1 + 'mod-' + year + '.tif')
x_grid, y_grid = np.meshgrid(mcd['x'], mcd['y'])

# Import ISMIP 1 km elevation data
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
elev = ismip_1km['SRF'].rio.write_crs('EPSG:3413')
elev_match = elev.rio.reproject_match(mcd)
mask = ismip_1km['GIMP'].rio.write_crs('EPSG:3413').astype(int)
mask_match = mask.rio.reproject_match(mcd)

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Define which MODIS data to use
modis = mcd

# Mask
modis = xr.where(mask_match, modis, np.nan)
sw = xr.where(mask_match, sw, np.nan)

#%%

# Define grid
grid = sw[:,:,0].values

# Define a 10x10 kernel for the sliding window
kernel = np.ones((10, 10), dtype=int)

# Apply 2D convolution to sum values in each 5x5 window
window_sum = convolve2d(grid, kernel, mode='same')

# Find locations where the sum of the 10x10 block is 0
zero_blocks = (window_sum == 0)

# Create a mask for NaN values
nan_mask = np.isnan(modis[:,:,0].values)

# Define a 5x5 kernel for the sliding window
kernel = np.ones((5, 5), dtype=int)

# Apply convolution to count NaNs in each 5x5 block
nan_count = convolve2d(nan_mask, kernel, mode="same")

# Find locations where the 5x5 block contains no NaNs (count == 0)
no_nan_blocks = (nan_count == 0)

#%%

# Define window size from 1 to 250 km
windows = np.array((1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90,
                    100))
kernel_size = (windows[:-1] + windows[1:])**2

final_data_list = []

for w in range(len(windows[:-1])):
    print('Window size is %.0f pixels' % ((windows[w] +  windows[w+1])**2))


    data_list = []
    
    for e in range(len(elevations) - 1):
        #print('Processing %.2f to %.2f m' %(bin_edges[i], bin_edges[i+1]))
        
        # Get first elevation band
        band = (elev_match > elevations[e]) & (elev_match < elevations[e+1]) &\
            (no_nan_blocks == True) & (mask_match.values == 1)
    
        # Mask
        masked_sw = xr.where(band == 1, zero_blocks, np.nan)
        
        # Find cells in band with 0% water
        x, y = np.where(masked_sw == 1)
        
        if x.shape[0] == 0:
            pass
        else:
        
            # Select number of samples
            np.random.seed(42)
            samples = np.random.choice(x.shape[0], size=sample_size, replace=False)
            #print('There are %.0f samples' % binned_counts[i])
            
            data_list1 = []
            
            for n in samples:
                
                # Get row/col
                coords = x[n], y[n]
                
                # Define window
                x0 = np.maximum(coords[0]-windows[w], 0)
                x1 = coords[0]+windows[w+1]
                y0 = np.maximum(coords[1]-windows[w], 0)
                y1 = coords[1]+windows[w+1]
                        
                # Identify nearest neighbors
                neighbors = np.array(sw[:,:,0][x0:x1,y0:y1]).flatten()
                albedo = np.array(modis[:,:,0][x0:x1,y0:y1]).flatten()
                
                # Find grid cells with more than 0% water
                indices = np.where(neighbors > 0)

                # Remove elements with more than 0% water
                albedo_filtered = np.where(neighbors > 0, np.nan, albedo)
                albedo_non_water = albedo_filtered[np.isfinite(albedo_filtered)]
                                
                # Append to list
                data_list1.append((y_grid[coords],x_grid[coords], np.nanstd(albedo_non_water),
                                  np.nanmean(albedo_non_water), elev_match[coords].values,
                                  albedo_non_water.shape[0]))
        
        # Append
        data_list.extend(data_list1)
    
    # Append
    final_data_list.append(data_list)
            

stats_list = []

for i in range(len(final_data_list)):
    # Make a DataFrame
    df = pd.DataFrame(final_data_list[i])
    df.columns = ['y', 'x', 'std', 'mean_albedo', 'elevation', 'count']
    df['bin'] = pd.cut(df['elevation'], bins=elevations, right=False)
    group = df.groupby('bin', observed=True)['std'].mean()
    stats_list.append(group.values)

#%%
# Make a DataFrame
df = pd.DataFrame(stats_list)

# Export
df.to_csv(path1 + 'scales/non.csv')


#%%
df = pd.DataFrame(final_data_list)
df.columns = ['y', 'x', 'std', 'mean_albedo', 'elevation']
df['elevation'] = df['elevation'].astype(float)
df.to_csv(path1 + 'non-water-albedo-effect-' + year + '.csv')

#%%



















