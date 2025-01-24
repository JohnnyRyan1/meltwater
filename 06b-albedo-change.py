#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Local albedo variability explained by surface water.

"""

# Import packages
import rioxarray as rio
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Define year
year = str(2019)

# Read raster data
sw = rio.open_rasterio(path1 + 'zhang/surface_water_mask_' + year +'_500m.tif')

# Replace NaN values with zeros
sw = sw.fillna(0)

mcd = rio.open_rasterio(path1 + 'mcd-' + year + '.tif')
mod = rio.open_rasterio(path1 + 'mod-' + year + '.tif')
x_grid, y_grid = np.meshgrid(mcd['x'], mcd['y'])

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
elev = ismip_1km['SRF'].rio.write_crs('EPSG:3413')
elev_match = elev.rio.reproject_match(mcd)
mask = ismip_1km['GIMP'].rio.write_crs('EPSG:3413').astype(int)
mask_match = mask.rio.reproject_match(mcd)

# Define which MODIS data to use
modis = mcd

# Mask
modis = xr.where(mask_match, modis, np.nan)
sw = xr.where(mask_match, sw, np.nan)


#%%

# Define window size of 3 km
window = 3

# Find cells with >50% water
x, y = np.where(sw[:,:,0] > 0.5)

data_list = []

for n in range(x.shape[0]):
    
    # Get row/col
    coords = x[n], y[n]
    
    # Define window
    x0 = np.maximum(coords[0]-window, 0)
    x1 = coords[0]+window+1
    y0 = np.maximum(coords[1]-window, 0)
    y1 = coords[1]+window+1
            
    # Identify nearest neighbors
    neighbors = np.array(sw[:,:,0][x0:x1,y0:y1]).flatten()
    albedo = np.array(modis[:,:,0][x0:x1,y0:y1]).flatten()
    
    # Only compute linear regression if number of albedo values is greater than 9
    if np.isfinite(albedo).sum() < 9:
        pass
    else:
        # Find grid cells where surface water is NaN
        indices = np.where(np.isnan(neighbors))

        # Remove elements with NaN surface water
        albedo_filtered = np.where(np.isnan(neighbors), np.nan, albedo)
        
        # Remove elements with NaN surface water AND albedo
        albedo_non_nan = albedo_filtered[np.isfinite(albedo_filtered)]
        water_non_nan = neighbors[np.isfinite(albedo_filtered)]

        # Linear regression
        slope, intercept, r, p, se = stats.linregress(water_non_nan, albedo_non_nan)
        
        # Append to list
        data_list.append((y_grid[coords],x_grid[coords], np.nanstd(albedo_non_nan),
                          np.nanmean(albedo_non_nan), intercept, slope, r, p, 
                          elev_match[coords].values, albedo_non_nan.shape[0]))


    
#%%

# Make a DataFrame
df = pd.DataFrame(data_list)
df.columns = ['y', 'x', 'std', 'mean_albedo', 'intercept', 'slope', 
              'r', 'p', 'elevation', 'grid_cells']
df['elevation'] = df['elevation'].astype(float)

# Filter
df = df[df['slope'] < 0]
df = df[df['p'] < 0.001]

df.to_csv(path1 + 'albedo-effects-' + year + '.csv')


#%%






