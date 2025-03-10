#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Meltwater runoff vs. energy absorbed by surface.

"""


# Import packages
import xarray as xr
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import matplotlib.ticker as ticker


#%%

# Define path
path1 = '/Users/jr555/Documents/research/hydrology/mar/'
path2 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/'
    
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
sw_energy_by_year = []
for file in mar_files:
    mar = xr.open_dataset(file)
    
    # Sum energy fluxes
    mar_sw = (mar['SWD']*(1-mar['AL2'][:,0,:,:]))[152:244,:,:].where(msk)
    swnet_energy = np.nansum(((np.nansum(mar_sw, axis=0) * 1000000 * 100) * 86400) / 1e15)
    
    mar_lw = mar['LWD'][152:244,:,:].where(msk)
    lw_energy = np.nansum(((np.nansum(mar_lw, axis=0) * 1000000 * 100) * 86400) / 1e15)
    
    mar_shf = mar['SHF'][152:244,:,:].where(msk)
    shf_energy = np.nansum(((np.nansum(mar_shf, axis=0) * 1000000 * 100) * 86400) / 1e15)
    
    mar_lhf = mar['LHF'][152:244,:,:].where(msk)
    lhf_energy = np.nansum(((np.nansum(mar_lhf, axis=0) * 1000000 * 100) * 86400) / 1e15)
    
    energy_by_year.append(swnet_energy+lw_energy+shf_energy+lhf_energy)
    sw_energy_by_year.append(swnet_energy)
    
    runoff_mm = np.nansum(mar['RU'][152:244,0,:,:].sum(axis=0).where(msk))
    runoff_gt = (runoff_mm * 0.001) * (100 * 1000000) / 1000000000
    melt_mm = np.nansum(mar['ME'][152:244,0,:,:].sum(axis=0).where(msk))
    melt_gt = (melt_mm * 0.001) * (100 * 1000000) / 1000000000
    runoff_by_year.append(runoff_gt)
    melt_by_year.append(melt_gt)

#%%

energy = sw_energy_by_year

# Compute best-fit line
coeffs = np.polyfit(energy, runoff_by_year, 1)
poly_eq = np.poly1d(coeffs)
r_value, _ = pearsonr(energy, runoff_by_year)
best_fit_line = poly_eq(energy)


# Calculate the residuals
residuals = runoff_by_year - best_fit_line

# Calculate the standard deviation of the residuals
std_residuals = np.std(residuals)

# Calculate the 90% confidence intervals
confidence_interval = 1.645 * std_residuals

# Generate x values for the confidence interval curves
x_values = np.linspace(min(energy)-50000, max(energy)+50000, 500)

# Generate y values for the confidence interval curves
y_values_upper = poly_eq(x_values) + confidence_interval * np.sqrt(1 + (x_values - np.mean(sw_energy_by_year))**2 / np.sum((sw_energy_by_year - np.mean(sw_energy_by_year))**2))
y_values_lower = poly_eq(x_values) - confidence_interval * np.sqrt(1 + (x_values - np.mean(sw_energy_by_year))**2 / np.sum((sw_energy_by_year - np.mean(sw_energy_by_year))**2))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot the area-frequency plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(5,5), 
                                    layout='constrained')

ax1.scatter(energy, runoff_by_year, label='', color=c1, lw=2, zorder=2)

# Plot best-fit line
ax1.plot(sorted(x_values), poly_eq(sorted(x_values)), color='black', 
        linestyle='dashed', linewidth=1)

# Plot curved confidence intervals
ax1.fill_between(x_values, y_values_lower, y_values_upper, color='gray', 
                 alpha=0.2, zorder=1)
ax1.set_ylabel("Summer runoff (Gt)", fontsize=12)
ax1.set_xlabel("Net summer shortwave radiation (PJ)", fontsize=12)

ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_xlim(720000, 850000)
ax1.set_ylim(75, 520)

# Use ScalarFormatter to format the x-axis
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-3, 3))
ax1.xaxis.set_major_formatter(formatter)

plt.savefig(path2 + 'sfig-4-energy-vs-runoff.png', dpi=300)

#%%






