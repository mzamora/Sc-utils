#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 11:38:24 2018

@author: elynn
"""
import sys
sys.path.append('/mnt/lab_45d1/users/elynn/sdge_operational_post_processing/')
import WRF_methods as wrf_m
import pygrib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.0)
filename = '/home/elynn/Documents/EDMF_downdraft/Validation/RAP_initial/rap.t09z.awp130bgrbf00.grib2'
f = pygrib.open(filename)
#for g in f:
#  if g.name=='Temperature':
#      print g.typeOfLevel#, g.level, g.name, g.validDate, g.analDate, g.forecastTime
  
lats, lons = g.latlons()
temperature = f.select(name='Temperature',typeOfLevel='hybrid')
qt = f.select(name='Specific humidity',typeOfLevel='hybrid')
ql = f.select(name='Cloud mixing ratio',typeOfLevel='hybrid')
P = f.select(name='Pressure',typeOfLevel='hybrid')
NKX_loc = wrf_m.find_nearest(lats,lons,32.85,-117.12)
NKX_temp = []
NKX_ql = []
NKX_qt = []
NKX_P = []
for i in range(len(temperature)):
    NKX_temp.append(temperature[i].values[NKX_loc[0][0],NKX_loc[1][0]])
    NKX_ql.append(ql[i].values[NKX_loc[0][0],NKX_loc[1][0]])
    NKX_qt.append(qt[i].values[NKX_loc[0][0],NKX_loc[1][0]])
    NKX_P.append(P[i].values[NKX_loc[0][0],NKX_loc[1][0]])
NKX_temp, NKX_ql, NKX_qt, NKX_P = np.array(NKX_temp), np.array(NKX_ql), np.array(NKX_qt), np.array(NKX_P)
cpAir = 1005.
R = 287. 
P0 = 100000.
NKX_theta =  NKX_temp * (P0/NKX_P)**(R/cpAir)
NKX_thetaV = NKX_theta*(1+0.61*NKX_qt-NKX_ql)
plt.figure(figsize=(7,5))
plt.subplot(1,3,1)
plt.plot(NKX_thetaV,np.arange(50),'-o')
plt.ylim([0,15])
plt.xlim([287,305])
plt.xlabel(r'$\theta_V$ [K]')
plt.ylabel('Hybrid level')
plt.subplot(1,3,2)
plt.plot(NKX_qt,np.arange(50),'-o')
plt.ylim([0,15])
plt.xlabel(r'$q_T$ [kg/kg]')
plt.subplot(1,3,3)
plt.plot(NKX_ql,np.arange(50),'-o')
plt.ylim([0,15])
plt.xlabel(r'$q_l$ [kg/kg]')
plt.tight_layout()
plt.savefig('/home/elynn/Documents/EDMF_downdraft/Validation/RAP_initial/RAP_2017052509Z.png',\
            dpi=200,bbox_inches='tight')