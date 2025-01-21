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
year = str(2018)

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'
path2 = '/Users/jr555/Documents/research/hydrology/'

# Read raster data
sw = rio.open_rasterio(path1 + 'zhang/surface_water_mask_' + year + '_500m.tif')

# Read MODIS product
mcd = rio.open_rasterio(path2 + 'modis/mcd-' + year + '.tif')
mod = rio.open_rasterio(path2 + 'modis/mod-' + year + '.tif')
x_grid, y_grid = np.meshgrid(mcd['x'], mcd['y'])

# Import ISMIP 1 km elevation data
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
elev = ismip_1km['SRF'].rio.write_crs('EPSG:3413')
elev_match = elev.rio.reproject_match(mcd)

# Import ISMIP 1 km ice mask
gimp_1km = rio.open_rasterio(path1 + '1km-ISMIP6-GIMP.tif')
gimp_1km = gimp_1km.rio.write_crs('EPSG:3413')
gimp_1km_match = gimp_1km.rio.reproject_match(mcd)

# Import coefficients
coeffs = pd.read_csv(path1 + 'albedo-effects-' + year + '.csv')

# Define sample distribution
sample = np.random.normal(coeffs['elevation'].mean(), coeffs['elevation'].std(), 
                          coeffs.shape[0])
bins = np.linspace(coeffs['elevation'].min(), coeffs['elevation'].max(), 20 + 1)
binned_counts, bin_edges = np.histogram(coeffs['elevation'], bins=bins)

#%%

modis = mcd

#%%

# Define grid
grid = sw[0,:,:].values

# Set values outside the ice sheet to NaN
grid[~gimp_1km_match[0,:,:].values.astype(bool)] = np.nan

# Define a 5x5 kernel for the sliding window
kernel = np.ones((10, 10), dtype=int)

# Apply 2D convolution to sum values in each 5x5 window
window_sum = convolve2d(grid, kernel, mode='same')

# Find locations where the sum of the 5x5 block is 0
zero_blocks = (window_sum == 0)

#%%

# Create a mask for NaN values
nan_mask = np.isnan(modis[0,:,:].values)

# Define a 5x5 kernel for the sliding window
kernel = np.ones((5, 5), dtype=int)

# Apply convolution to count NaNs in each 5x5 block
nan_count = convolve2d(nan_mask, kernel, mode="same")

# Find locations where the 5x5 block contains no NaNs (count == 0)
no_nan_blocks = (nan_count == 0)

#%%
final_data_list = []

for i in range(len(bins) - 1):
    print('Processing %.2f mto %.2f m' %(bin_edges[i], bin_edges[i+1]))
    
    # Get first elevation band
    band = (elev_match > bin_edges[i]) & (elev_match < bin_edges[i+1]) &\
        (no_nan_blocks == True) & (gimp_1km_match[0,:,:].values == 1)

    # Mask
    masked_sw = xr.where(band == 1, zero_blocks, np.nan)
    
    # Find cells in band with 0% water
    x, y = np.where(masked_sw == 1)
    
    if x.shape[0] == 0:
        pass
    else:
    
        # Select number of samples
        np.random.seed(42)
        samples = np.random.choice(x.shape[0], size=binned_counts[i], replace=False)
        print('There are %.0f samples' % binned_counts[i])
        
        data_list = []
        
        for n in samples:
            
            # Get row/col
            coords = x[n], y[n]
            
            # Identify nearest neighbors (5 x 5)
            neighbors = np.array(sw[0,:,:][coords[0]-2:coords[0]+3,
                                                coords[1]-2:coords[1]+3]).flatten()
            
            # Confirm that there is no surface water
            if np.sum(neighbors) > 0:
                print("Error...")
            else:
                albedo = np.array(modis[0,:,:][coords[0]-2:coords[0]+3,
                                             coords[1]-2:coords[1]+3]).flatten()
                
                # Append to list
                data_list.append((y_grid[coords],x_grid[coords], np.nanstd(albedo),
                                  np.nanmean(albedo), elev_match[coords].values))
    
    # Append
    final_data_list.extend(data_list)
            
#%%
df = pd.DataFrame(final_data_list)
df.columns = ['y', 'x', 'std', 'mean_albedo', 'elevation']
df['elevation'] = df['elevation'].astype(float)
df.to_csv(path1 + 'non-water-albedo-effect-' + year + '.csv')

#%%



















