#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Meltwater runoff vs. energy absorbed by surface.

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
path1 = '/Users/jr555/Documents/research/hydrology/mar/'

# Define files
mar_files = sorted(glob.glob(path1 + 'MARv3.12.1-10km-daily-ERA5-*.nc'))

# Derive MAR surface heights
mar = xr.open_dataset(mar_files[0])
sh = mar['SH'].values
msk = np.isfinite(mar['MSK'].where(mar['MSK'] > 90))

#%%

# Read first file
energy_by_year = []
runoff_by_year = []
melt_by_year = []
for file in mar_files:
    mar = xr.open_dataset(file)
    watts = mar['LWD'] + (mar['SWD']*(1-mar['AL2'][:,0,:,:])) + mar['SHF'] + mar['LHF']
    watts_sum = np.nansum(watts.where(msk)[152:244,:,:].sum(axis=0))
    
    watts_sum * 1000000 * 100 * 86400 / 1e15
    
    runoff_mm = np.nansum(mar['RU'][152:244,0,:,:].sum(axis=0).where(msk))
    runoff_gt = (runoff_mm * 0.001) * (100 * 1000000) / 1000000000
    melt_mm = np.nansum(mar['ME'][152:244,0,:,:].sum(axis=0).where(msk))
    melt_gt = (melt_mm * 0.001) * (100 * 1000000) / 1000000000
    #energy_by_year.append()
    runoff_by_year.append(runoff_gt)
    melt_by_year.append(melt_gt)

#%%

mar_sw = (mar['SWD']*(1-mar['AL2'][:,0,:,:]))[152:244,:,:].where(msk)
swnet_energy = np.nansum(((np.nanmean(mar_sw, axis=0) * 1000000 * 100) * 86400) / 1e15)

mar_lw = mar['LWD'].where(msk)
lw_energy = np.nansum(((np.nanmean(mar_lw, axis=0) * 1000000 * 100) * 86400) / 1e15)

mar_shf = mar['SHF'].where(msk)
shf_energy = np.nansum(((np.nanmean(mar_shf, axis=0) * 1000000 * 100) * 86400) / 1e15)

mar_lhf = mar['LHF'].where(msk)
lhf_energy = np.nansum(((np.nanmean(mar_lhf, axis=0) * 1000000 * 100) * 86400) / 1e15)


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
ax1.plot(times, np.cumsum(runoff[1][:-1]), label='2016', color=c2)
ax1.plot(times, np.cumsum(runoff[2]), label='2018', color=c4, lw=2)
ax1.plot(times, np.cumsum(runoff[3]), label='2019', color=c1, lw=2)

ax1.vlines(x=pd.to_datetime(drone_2015), ymin=0, ymax=np.cumsum(runoff[0])[idx1][0], 
            ls='dashed', lw=1.5, color=c3)
ax1.hlines(y=np.cumsum(runoff[0])[idx1][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(drone_2015), ls='dashed', lw=1.5, color=c3)

ax1.vlines(x=pd.to_datetime(z_2018), ymin=0, ymax=np.cumsum(runoff[2])[idx3][0], 
            ls='dashed', lw=1.5, color=c4)
ax1.hlines(y=np.cumsum(runoff[2])[idx3][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(z_2018), ls='dashed', lw=1.5, color=c4)

ax1.vlines(x=pd.to_datetime(z_2019), ymin=0, ymax=np.cumsum(runoff[3])[idx4][0], 
            ls='dashed', lw=1.5, color=c1)
ax1.hlines(y=np.cumsum(runoff[3])[idx4][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(z_2019), ls='dashed', lw=1.5, color=c1)

ax1.vlines(x=pd.to_datetime(drone_2016), ymin=0, ymax=np.cumsum(runoff[1])[idx2][0], 
            ls='dashed', lw=1.5, color=c2)
ax1.hlines(y=np.cumsum(runoff[1])[idx2][0], xmin=pd.to_datetime(start_date), 
           xmax=pd.to_datetime(drone_2016), ls='dashed', lw=1.5, color=c2)

ax1.set_ylabel("Cumulative runoff (mm w.e.)", fontsize=12)
ax1.legend(fontsize=12, loc=4)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
ax1.set_ylim(0, 1850)

#%%






