#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make a plot showing accuracy scores vs. NDWI thresholds.

"""

# Import packages
import pandas as pd
import matplotlib.pyplot as plt

#%%

# Read
df_2015 = pd.read_csv('/Users/jr555/Documents/research/hydrology/ndwi-thresholds-2015.csv')
df_2016 = pd.read_csv('/Users/jr555/Documents/research/hydrology/ndwi-thresholds-2016.csv')

#%%

# Plot
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12,4), 
                                    layout='constrained')

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(df_2015['threshold'], df_2015['precision'], zorder=2, color=c3, lw=2, marker='o', linestyle='-',
         label='Precision')
ax1.plot(df_2015['threshold'], df_2015['accuracy'], zorder=2, color=c2, lw=2, marker='o', linestyle='-',
         label='Accuracy')
ax1.legend(fontsize=12)
ax1.set_xlabel('NDWI threshold', fontsize=12)
ax1.set_ylabel('Accuracy', fontsize=12)
ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=12)


ax2.plot(df_2016['threshold'], df_2016['precision'], zorder=2, color=c3, lw=2, marker='o', linestyle='-',
         label='Precision')
ax2.plot(df_2016['threshold'], df_2016['accuracy'], zorder=2, color=c2, lw=2, marker='o', linestyle='-',
         label='Accuracy')
ax2.legend(fontsize=12)
ax2.set_xlabel('NDWI threshold', fontsize=12)
ax2.set_ylabel('Accuracy', fontsize=12)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=12)

#%%





























