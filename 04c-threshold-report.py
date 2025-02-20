#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots showing choice of NDWI threshold.

"""

# Import packages
import pandas as pd
import matplotlib.pyplot as plt

# Import data
russell_df = pd.read_csv('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/ndwi-thresholds-russell.csv')
wq7_df = pd.read_csv('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/ndwi-thresholds-2015.csv')

#%%
# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(7,4), 
                                    layout='constrained')

ax1.plot(wq7_df['threshold'], wq7_df['recall'], lw=2, zorder=0,
         label='Recall', color=c1)
ax1.plot(wq7_df['threshold'], wq7_df['accuracy'], lw=2, zorder=0, 
         label='Accuracy', color=c2)
ax1.plot(wq7_df['threshold'], wq7_df['precision'], lw=2, zorder=0,
         label='Precision', color=c3)

ax1.set_xlabel("NDWI threshold", fontsize=12)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_ylabel("Score", fontsize=12)
ax1.set_xlim(0, 0.3)
ax1.set_ylim(0.4, 1.03)
ax1.legend(loc=3, fontsize=12)

#%%

# Plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(7,4), 
                                    layout='constrained')

ax1.plot(russell_df['threshold'], russell_df['recall'], lw=2, zorder=0,
         label='Recall', color=c1)
ax1.plot(russell_df['threshold'], russell_df['accuracy'], lw=2, zorder=0, 
         label='Accuracy', color=c2)
ax1.plot(russell_df['threshold'], russell_df['precision'], lw=2, zorder=0,
         label='Precision', color=c3)

ax1.set_xlabel("NDWI threshold", fontsize=12)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_ylabel("Score", fontsize=12)
ax1.set_xlim(0, 0.3)
ax1.set_ylim(0.4, 1.03)
ax1.legend(loc=3, fontsize=12)


#%%

