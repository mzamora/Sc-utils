#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 15:17:41 2018

@author: elynn
"""
from netCDF4 import Dataset
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks",font_scale=1.1)
def draw_z_grid_background(ax,zm,xmin,xmax):
    for z in range(len(zm)):
        ax.plot([xmin,xmax],[zm[z],zm[z]],color='gray',alpha=0.7,lw=0.5)
radtyp4 = Dataset('/mnt/lab_45d1/database/Sc_group/uclales_output/archive/2xPBL_or_1km/NKX_20170511/hr_4_18_60min/NKX_20170511.ps.nc')
radtyp2 = Dataset('/mnt/lab_45d1/database/Sc_group/uclales_output/NKX_20170511/hr_4_18_60min/NKX_20170511.ps.nc')
thetaL_4 = radtyp4['lflxu'][:,:]-radtyp4['lflxd'][:,:]
thetaL_2 = radtyp2['rflx']#radtyp2['lflxu'][:,:]-radtyp2['lflxd'][:,:]
time = radtyp2['time'][:]
z_2 = radtyp2['zm'][:]  
z_4 = radtyp4['zm'][:]
t_index = 1
plt.figure(figsize=(5,8))
ax1 = plt.subplot(111)
colors2 = sns.cubehelix_palette(8, start=2, rot=0.05, dark=0, light=.95, reverse=True)
draw_z_grid_background(ax1,z_4,thetaL_2[55,:].min()-1,thetaL_4[55,:].max()+1)
plt.plot(thetaL_2[t_index,:],z_2,'--',color='blue')
plt.plot(thetaL_4[t_index,:],z_4,'--',color='red')
plt.plot(thetaL_2[55,:],z_2,color='blue',label='radtype 2')
plt.plot(thetaL_4[55,:],z_4,color='red',label='radtype 4')
plt.legend(loc=0)
#plt.xlim([thetaL_2[t_index,:].min()-1,310.])#thetaL_2[t_index,:].max()+1])
plt.ylim([0,2000])
plt.xlabel(r'Net LW radiative flux $[W/m^2]$')
#plt.xlabel(r'$\theta_l [K]$')
plt.ylabel(r'z [m]')
#plt.savefig('/mnt/lab_45d1/database/Sc_group/LES_analysis/radtyp2_vs_4/NKX_20170511_net_LW_hr6_and_hr20.png',dpi=200,bbox_inches='tight')