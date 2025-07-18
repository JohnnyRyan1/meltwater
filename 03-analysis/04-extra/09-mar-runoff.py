#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cumulative runoff at WQ7 + Russell + Isunguata Sermia from MAR

Zhang et al. (2022):
    Period for 2018 is Jul 25 - Jul 30
    Period for 2019 is Jul 29 - Aug 5

"""


# Import packages
import xarray as xr
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


#%%

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/mar/'
path2 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/'

# Define files
mar_files = sorted(glob.glob(path1 + 'MARv3.12.1-10km-daily-ERA5-*.nc'))

# Define dates
drone_2015 = "2019-07-21"
drone_2016 = "2019-07-05"
z_2018 = "2019-07-27"
z_2019 = "2019-08-02"

#%%

# Read first file
mar = xr.open_dataset(mar_files[2])

# Find grid cells for WQ7
target_lat = 67.05
target_lon = -48.9

# Calculate the distance to the target point
distances = np.sqrt((mar['LAT'].values - target_lat)**2 + (mar['LON'].values - target_lon)**2)

# Find the index of the minimum distance
min_dist_idx = np.unravel_index(np.argmin(distances), mar['LAT'].values.shape)

#%%

times = pd.to_datetime(mar['TIME'].values).normalize()
idx1 = np.where(times == pd.Timestamp(drone_2015))[0]
idx2 = np.where(times == pd.Timestamp(drone_2016))[0]
idx3 = np.where(times == pd.Timestamp(z_2018))[0]
idx4 = np.where(times == pd.Timestamp(z_2019))[0]

runoff = []
for file in mar_files:
    
    # Read first file
    mar = xr.open_dataset(file)

    # Define runoff
    runoff.append(mar['RU'][:,0,min_dist_idx[0], min_dist_idx[1]].values)

#%%

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

start_date = "2019-06-01"
end_date = "2019-10-01"

# Plot the area-frequency plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8,4), 
                                    layout='constrained')

ax1.plot(times, np.cumsum(runoff[0]), label='2015', color=c3, lw=2)
#ax1.plot(times, np.cumsum(runoff[1]), label='2016', color=c2)
ax1.plot(times, np.cumsum(runoff[1]), label='2018', color=c4, lw=2)
ax1.plot(times, np.cumsum(runoff[2]), label='2019', color=c1, lw=2)

ax1.vlines(x=pd.to_datetime(drone_2015), ymin=0, ymax=np.cumsum(runoff[0])[idx1][0], 
            ls='dashed', lw=1.5, color=c3)
ax1.hlines(y=np.cumsum(runoff[0])[idx1][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(drone_2015), ls='dashed', lw=1.5, color=c3)

ax1.vlines(x=pd.to_datetime(z_2018), ymin=0, ymax=np.cumsum(runoff[1])[idx3][0], 
            ls='dashed', lw=1.5, color=c4)
ax1.hlines(y=np.cumsum(runoff[1])[idx3][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(z_2018), ls='dashed', lw=1.5, color=c4)

ax1.vlines(x=pd.to_datetime(z_2019), ymin=0, ymax=np.cumsum(runoff[2])[idx4][0], 
            ls='dashed', lw=1.5, color=c1)
ax1.hlines(y=np.cumsum(runoff[2])[idx4][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(z_2019), ls='dashed', lw=1.5, color=c1)

#ax1.vlines(x=pd.to_datetime(drone_2016), ymin=0, ymax=np.cumsum(runoff[1])[idx2][0], 
#            ls='dashed', lw=1.5, color=c2)
#ax1.hlines(y=np.cumsum(runoff[1])[idx2][0], xmin=pd.to_datetime(start_date), 
#           xmax=pd.to_datetime(drone_2016), ls='dashed', lw=1.5, color=c2)

ax1.set_ylabel("Cumulative runoff (mm w.e.)", fontsize=12)
ax1.legend(fontsize=12, loc=4)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
ax1.set_ylim(0, 1850)

plt.savefig(path2 + 'sfig-3-mar-runoff.png', dpi=300)

#%%






