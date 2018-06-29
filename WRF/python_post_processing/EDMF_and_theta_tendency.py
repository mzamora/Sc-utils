#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 14:16:03 2018

@author: elynn
"""

import sys
sys.path.append('/mnt/lab_45d1/database/Sc_group/github/Sc-utils/WRF/python_post_processing/')
import WRF_methods as wrf_m
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import datetime
from netCDF4 import Dataset
from wrf import getvar, destagger
from mpl_toolkits.basemap import Basemap
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.0)
rho = 1.28
Lv = 2.5e6
cpAir = 1005. #specific heat of dry air [J/kgK]
current_palette = sns.color_palette('deep')
'''
MYNN-RAD real simulation
'''
#WRF_config = ['MYNN-RAD']
#wrf_dirs = ['/mnt/lab_45d1/database/Sc_group/WRF_39_RAP/MYNN_RAD/wrfout_withSpinup/']
#hh = 12
#current_month = 5
#current_day = 25
#wrf_dir = wrf_dirs[0]
#time = pd.date_range('2017-'+str(current_month)+'-'+str(current_day)+' 10:15','2017-'+str(current_month)+'-'+str(current_day+1)+' 6:00',freq='15min')
#wrf_dir = wrf_dir+time[0].strftime('%Y%m%d')+'/'
#f = Dataset(wrf_dir+'wrfout_d01_2017'+time[0].strftime('%m%d')+'_10_15_00')
#lat = f['XLAT'][0,:,:]
#lon = f['XLONG'][0,:,:]
#NKX_loc = wrf_m.find_nearest(lat,lon,32.905680, -118.141179)#32.85,-117.12)
#time_pd = pd.DataFrame(np.arange(len(time)),index=time,columns=['index'])
#t_index = 7 #fixed at 12Z
#z = getvar(f,'z',timeidx=t_index)[:,NKX_loc[0][0],NKX_loc[1][0]]
#'''EDMF - heat flux'''
#thetaL = wrf_m.get_thetaL_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)
#MF_wpthetap = f['DBG5'][t_index,:,NKX_loc[0][0],NKX_loc[1][0]] #MF, already in W/m^2
#MF_wpthetap = MF_wpthetap/rho/cpAir #convert back to m*K/s
#ED_wpthetap, k1 = wrf_m.get_eddy_diffusivity_flux(f,z,t_index,NKX_loc[0][0],NKX_loc[1][0],thetaL,True)
#dEDMF_dz = wrf_m.dArray_dz_2nd_central(MF_wpthetap,z) + wrf_m.dArray_dz_2nd_central(ED_wpthetap,z)
#output_dtheta_dt = f['RTHBLTEN'][t_index+1,:,NKX_loc[0][0],NKX_loc[1][0]]
#plt.plot(output_dtheta_dt,z,'-o',color='k',label='Actual')
#plt.plot(-dEDMF_dz,z,'-o',color='b',label='Calculated')
#plt.ylim([0,2000.])

'''
MYNN-RAD WRF ideal (ocean only)
'''
filename = '/mnt/lab_48tb1/users/hyang/database/wrfIdeal/rap13km/dycomsrf01/results_rapres_downdraft/HYDD_idealTest_v1_GB2001_confirm_min_outside_noUpdraft_zi_05a_1.5ms_edmf1_pbl5_sw4_lw4_icloud1_icloudpbl0_ml2/wrfout_d01_2015-07-10_08_00_00'
f = Dataset(filename)
time = pd.date_range('2017-7-10 8:00','2017-7-10 14:00',freq='15min')
lat = f['XLAT'][0,:,:]
lon = f['XLONG'][0,:,:]
loc = wrf_m.find_nearest(lat,lon,32.905680, -123.141179)
time_pd = pd.DataFrame(np.arange(len(time)),index=time,columns=['index'])
t_index = 30
z = getvar(f,'z',timeidx=t_index)[:,loc[0][0],loc[1][0]]
thetaL = wrf_m.get_thetaL_at_time_index(f,t_index,loc[0][0],loc[1][0],False)
MF_wpthetap = f['DBG5'][t_index,:,loc[0][0],loc[1][0]] #MF, already in W/m^2
MF_wpthetap = MF_wpthetap/rho/cpAir #convert back to m*K/s
ED_wpthetap, k1 = wrf_m.get_eddy_diffusivity_flux(f,z,t_index,loc[0][0],loc[1][0],thetaL,True)
dEDMF_dz = wrf_m.dArray_dz_2nd_central(MF_wpthetap,z) + wrf_m.dArray_dz_2nd_central(ED_wpthetap,z)
output_dtheta_dt = f['RTHBLTEN'][t_index+1,:,loc[0][0],loc[1][0]]
plt.plot(output_dtheta_dt,z,'-o',color='k',label='Actual')
plt.plot(-dEDMF_dz,z,'-o',color='b',label='Calculated')
plt.ylim([0,2000.])
