#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 12:40:30 2025

@author: jr555
"""

import pandas as pd
import matplotlib.pyplot as plt

path = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/fig1/'

ne = pd.read_csv(path + 'ne1.csv')
n = pd.read_csv(path + 'n1.csv')
sw = pd.read_csv(path + 'sw1.csv')

#%%

# Plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(10, 4),
                                             layout='constrained', sharex=True)

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'


ax1.plot(n['distance']/1000, n['albedo'], lw=1.5, color=c1, zorder=3, label='Albedo')
ax1.grid(ls='dashed', lw=1, zorder=1)
ax1.set_ylabel('Albedo', fontsize=12)    
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_ylim(0.21, 0.75)
ax1.set_xlim(0, 150)
ax1.legend(loc=3)

ax2.plot(ne['distance']/1000, ne['albedo'], lw=1.5, color=c1, zorder=2)
ax2.grid(ls='dashed', lw=1, zorder=1)
ax2.set_ylabel('Albedo', fontsize=12)    
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.set_ylim(0.21, 0.75)
ax2.set_xlim(0, 150)

ax3.plot(sw['distance']/1000, sw['albedo'], lw=1.5, color=c1, zorder=2)
ax3.grid(ls='dashed', lw=1, zorder=1)
ax3.set_ylabel('Albedo', fontsize=12)    
ax3.set_xlabel('Distance (km)', fontsize=12)  
ax3.tick_params(axis='both', which='major', labelsize=12)
ax3.set_ylim(0.21, 0.75)
ax3.set_xlim(0, 150)

ax4 = ax1.twinx()
ax4.fill_between(n['distance']/1000, n['water'], lw=1, color=c4, zorder=2,
                 alpha=0.5, label='Water')
ax4.tick_params(axis='both', which='major', labelsize=12)
ax4.set_ylim(0, 1)
ax4.legend(loc=3)

ax5 = ax2.twinx()
ax5.fill_between(ne['distance']/1000, ne['water'], lw=1, color=c4, zorder=2,
                 alpha=0.5, label='Water')
ax5.tick_params(axis='both', which='major', labelsize=12)
ax5.set_ylim(0, 1)

ax6 = ax3.twinx()
ax6.fill_between(sw['distance']/1000, sw['water'], lw=1, color=c4, zorder=2,
                 alpha=0.5, label='Water')
ax6.tick_params(axis='both', which='major', labelsize=12)
ax6.set_ylim(0, 1)

plt.savefig(path + 'all-transects.svg')

#%%






#%%



