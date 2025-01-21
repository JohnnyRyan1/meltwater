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
path2 = '/Users/jr555/Documents/research/hydrology/'

# Define year
year = str(2018)

# Read raster data
sw = rio.open_rasterio(path1 + 'zhang/surface_water_mask_' + year +'_500m.tif')

mcd = rio.open_rasterio(path2 + 'modis/mcd-' + year + '.tif')
mod = rio.open_rasterio(path2 + 'modis/mod-' + year + '.tif')
x_grid, y_grid = np.meshgrid(mcd['x'], mcd['y'])

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
elev = ismip_1km['SRF'].rio.write_crs('EPSG:3413')
elev_match = elev.rio.reproject_match(mcd)

# Define which MODIS data to use
modis = mcd

#%%

# Find cells with >50% water
x, y = np.where(sw[0,:,:] > 0.5)

data_list = []
albedo_list = []
water_list = []

for n in range(x.shape[0]):
    
    # Get row/col
    coords = x[n], y[n]
    
    # Identify nearest neighbors (5 x 5)
    neighbors = np.array(sw[0,:,:][coords[0]-2:coords[0]+3,
                                        coords[1]-2:coords[1]+3]).flatten()
    albedo = np.array(modis[0,:,:][coords[0]-2:coords[0]+3,
                                 coords[1]-2:coords[1]+3]).flatten()
    
    # Sort
    sorted_albedo = [x for _, x in sorted(zip(neighbors, albedo))]
    sorted_water = sorted(neighbors)
    
    if np.isnan(sorted_albedo).any() | np.isnan(sorted_water).any():
        pass
    else:
        # Subrtract max value so albedo change is negative
        albedo_normal = sorted_albedo - np.nanmax(sorted_albedo)
        albedo_list.append(albedo_normal)
        water_list.append(sorted_water)
        
        # Linear regression
        slope, intercept, r, p, se = stats.linregress(sorted_water, albedo_normal)
        
        # Append to list
        data_list.append((y_grid[coords],x_grid[coords], np.nanstd(sorted_albedo),
                          np.nanmean(sorted_albedo), intercept, slope, r, p, 
                          elev_match[coords].values))

#%%
# Make a DataFrame
df = pd.DataFrame(data_list)
df.columns = ['y', 'x', 'std', 'mean_albedo', 'intercept', 'slope', 'r', 'p', 'elevation']
df['elevation'] = df['elevation'].astype(float)

# Filter
df = df[df['slope'] < 0]
df = df[df['p'] < 0.001]

df.to_csv(path1 + 'albedo-effects-' + year + '.csv')


#%%






