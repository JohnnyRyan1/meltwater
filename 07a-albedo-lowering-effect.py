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
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'
path2 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/'

# Define year
year = str(2018)

# Read raster data
sw = rio.open_rasterio(path1 + 'zhang/surface_water_mask_' + year +'_500m.tif')

# Replace NaN values with zeros
sw = sw.fillna(0)

mcd = rio.open_rasterio(path1 + 'mcd-' + year + '.tif')
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

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Define window size from 1 to 100 km
windows = np.array((2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 25, 50, 100))

# Find cells with >50% water
x, y = np.where(sw[:,:,0] > 0.5)

#%%

# Define lists
final_data_list = []

for w in range(len(windows[:-1])):
    
    print('Window size is %.0f pixels' % ((windows[w] +  windows[w+1])**2))

    data_list = []

    for n in range(x.shape[0]):
        
        # Get row/col
        coords = x[n], y[n]
        
        # Define window
        x0 = np.maximum(coords[0]-windows[w], 0)
        x1 = coords[0]+windows[w]+1
        y0 = np.maximum(coords[1]-windows[w], 0)
        y1 = coords[1]+windows[w]+1
                
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
            data_list.append((y_grid[coords],x_grid[coords], 
                              np.nanmean(albedo_non_nan), intercept, slope, r, p, se,
                              elev_match[coords].values, albedo_non_nan.shape[0]))
    
    final_data_list.append(data_list)

# Export individual DataFrames for each length scale
for i in range(len(final_data_list)):
    # Make a DataFrame
    df = pd.DataFrame(final_data_list[i])
    df.columns = ['y', 'x', 'mean_albedo', 'intercept', 'slope', 
                  'r', 'p', 'se', 'elevation', 'grid_cells']
    # Filter
    df = df[df['slope'] < 0]
    
    # Export
    df.to_csv(path1 + 'effects/albedo-effects-' + year + '-' + str(windows[i]).zfill(2) + 'km.csv', index=False)
    

#%%
stats_list_slope = []
stats_list_r = []
stats_list_p = []
stats_list_se = []

for i in range(len(final_data_list)):
    # Make a DataFrame
    df = pd.DataFrame(final_data_list[i])
    df.columns = ['y', 'x', 'mean_albedo', 'intercept', 'slope', 
                  'r', 'p', 'se', 'elevation', 'grid_cells']
    # Filter
    df = df[df['slope'] < 0]
    
    df['bin'] = pd.cut(df['elevation'], bins=elevations, right=False)
    group_slope = df.groupby('bin', observed=True)['slope'].mean()
    stats_list_slope.append(group_slope.values)
    group_r = df.groupby('bin', observed=True)['r'].mean()
    stats_list_r.append(group_r.values)
    group_p = df.groupby('bin', observed=True)['p'].mean()
    stats_list_p.append(group_p.values)
    group_se = df.groupby('bin', observed=True)['se'].mean()
    stats_list_se.append(group_se.values)

# Make a DataFrame
df_slope = pd.DataFrame(stats_list_slope)
df_r = pd.DataFrame(stats_list_r)
df_p = pd.DataFrame(stats_list_p)
df_se = pd.DataFrame(stats_list_se)

# Export
df_slope.to_csv(path1 + 'effects/albedo-effects-slope-' + year + '.csv', index=False)
df_r.to_csv(path1 + 'effects/albedo-effects-r-' + year + '.csv', index=False)
df_p.to_csv(path1 + 'effects/albedo-effects-p-' + year + '.csv', index=False)
df_se.to_csv(path1 + 'effects/albedo-effects-se-' + year + '.csv', index=False)

#%%

# Sample 20 grid cells for plot at length scales of 5km

w=3

# Select number of samples
np.random.seed(42)
samples = np.random.choice(x.shape[0], size=50, replace=False)

water_list = []
albedo_list = []

for n in samples:
    
    # Get row/col
    coords = x[n], y[n]
            
    # Define window
    x0 = np.maximum(coords[0]-windows[w], 0)
    x1 = coords[0]+windows[w]+1
    y0 = np.maximum(coords[1]-windows[w], 0)
    y1 = coords[1]+windows[w]+1
            
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
        
        if r < -0.5:
            # Append to list
            water_list.append(water_non_nan)
            albedo_list.append(albedo_non_nan)
        else:
            pass
        
#%%

# Plot showing correlation coefficient with length scale
fig, axes = plt.subplots(nrows=5, ncols=4, 
                         figsize=(8, 8), layout='constrained',
                         sharex=True, sharey=True)

# Flatten the 2D array of axes to iterate easily
axes = axes.flatten()

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Loop through each subplot and plot data
for i, ax in enumerate(axes):
    ax.scatter(water_list[i], albedo_list[i], s=25)
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Compute best-fit line
    coeffs = np.polyfit(water_list[i], albedo_list[i], 1)
    poly_eq = np.poly1d(coeffs)
    r_value, _ = pearsonr(water_list[i], albedo_list[i])
    
    # Plot best-fit line
    ax.plot(sorted(water_list[i]), poly_eq(sorted(water_list[i])), color='black', 
            linestyle='dashed', linewidth=1)
    
    # Add r-value text to subplot
    ax.text(0.05, 0.05, f'r = {r_value:.2f}', transform=ax.transAxes, fontsize=11, 
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='gray'))
    ax.set_xlim(0, 1)

fig.supxlabel('Fraction of meltwater ponding', fontsize=12)
fig.supylabel('Albedo', fontsize=12)

plt.savefig(path2 + 'sfig-1-water-vs-albedo.png', dpi=300)


#%%





