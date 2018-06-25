#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 14:03:50 2018

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
WRF_config = ['MYNN-RAD']#,'MYNN-STEM','YSU','MYNN-RAD','MYNN-STEMRAD']
wrf_dirs = ['/mnt/lab_45d1/database/Sc_group/WRF_39_RAP/MYNN_RAD/wrfout_withSpinup/']
months = [5]#[5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,9,9,9]
days = [25]#[11,12,25,29,30,3,4,7,8,23,28,2,4,5,12,14,17,22,1,6,8,9,11,13,14,22,24,5,7,20]
hh = 12
current_month = 5
current_day = 25
wrf_dir = wrf_dirs[0]
time = pd.date_range('2017-'+str(current_month)+'-'+str(current_day)+' 10:15','2017-'+str(current_month)+'-'+str(current_day+1)+' 6:00',freq='15min')
wrf_dir = wrf_dir+time[0].strftime('%Y%m%d')+'/'
f = Dataset(wrf_dir+'wrfout_d01_2017'+time[0].strftime('%m%d')+'_10_15_00')
lat = f['XLAT'][0,:,:]
lon = f['XLONG'][0,:,:]
NKX_loc = wrf_m.find_nearest(lat,lon,32.85,-117.12)
time_pd = pd.DataFrame(np.arange(len(time)),index=time,columns=['index'])
t_index = time_pd.between_time('11:15','12:00').values[:,0]
z = getvar(f,'z',timeidx=t_index[0])[:,NKX_loc[0][0],NKX_loc[1][0]]

'''EDMF - heat flux'''
thetaL = wrf_m.get_thetaL_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)
MF_wpthetap = f['DBG5'][t_index,:,NKX_loc[0][0],NKX_loc[1][0]] #MF, already in W/m^2
MF_wpthetap = np.average(MF_wpthetap,axis=0)
ED_wpthetap, k1 = wrf_m.get_eddy_diffusivity_flux(f,z,t_index,NKX_loc[0][0],NKX_loc[1][0],thetaL,True)
ED_wpthetap, k1 = rho*cpAir*ED_wpthetap, k1
ql = wrf_m.get_qL_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)

'''EDMF - moisture flux'''
qt = wrf_m.get_qt_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)
MF_wpqvp = f['DBG22'][t_index,:,NKX_loc[0][0],NKX_loc[1][0]] #MF, already in W/m^2
MF_wpqvp = np.average(MF_wpqvp,axis=0)
MF_wpqlp = f['DBG17'][t_index,:,NKX_loc[0][0],NKX_loc[1][0]] #MF, already in W/m^2
MF_wpqlp = np.average(MF_wpqlp,axis=0)
ED_wpqtp, k2 = wrf_m.get_eddy_diffusivity_flux(f,z,t_index,NKX_loc[0][0],NKX_loc[1][0],qt,True)
ED_wpqtp, k2 = rho*Lv*ED_wpqtp, k2

###plotting begins here
fig = plt.figure(figsize=(10,6))
gs = GridSpec(1,5)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])  
ax3 = fig.add_subplot(gs[0,2])  
ax4 = fig.add_subplot(gs[0,3])  
ax5 = fig.add_subplot(gs[0,4])  

ax1.plot(np.average(thetaL,axis=0),z,'-o')
ax1.set_ylim([0,2000])
ax1.set_xlabel(r'$\theta_l$ [K]')
ax1.set_xlim([int(thetaL.min())-0.5,305])

ax2.plot(np.average(ql,axis=0)*1000.,z,'-o')
ax2.set_ylim([0,2000])
ax2.set_xlabel(r'$q_l$ [g/kg]')

ax3.plot(np.average(qt,axis=0)*1000.,z,'-o')
ax3.set_ylim([0,2000])
ax3.set_xlabel(r'$q_T$ [g/kg]')

ax4.plot(ED_wpthetap,z[0:-1],'-o',label='ED',color='saddlebrown')
ax4.plot(MF_wpthetap,z,'-o',label='MF',color='navy')
ax4.plot(MF_wpthetap[0:-1]+ED_wpthetap,z[0:-1],'-o',label='Total',color='gray')
ax4.set_ylim([0,2000])
ax4.set_xlabel(r'$\rho c_p \overline{w^\prime {\theta_L^\prime}}\/\/[W/m^2]$')

ax5.plot(ED_wpqtp,z[0:-1],'-o',label='ED',color='saddlebrown')
ax5.plot(MF_wpqlp+MF_wpqvp,z,'-o',label='MF',color='navy')
ax5.plot(MF_wpqlp[0:-1]+MF_wpqvp[0:-1]+ED_wpqtp,z[0:-1],'-o',label='Total',color='gray')
ax5.set_ylim([0,2000])
ax5.set_xlabel(r'$\rho L \overline{w^\prime {q_t^\prime}}\/\/[W/m^2]$')
ax5.legend(loc=0)
for i in range(1,5):
    plt.setp(fig.axes[i].get_yticklabels(), visible=False)
plt.tight_layout()
f.close()
###Now read regular MYNN
wrf_dir = '/mnt/lab_45d1/database/Sc_group/WRF_39_RAP/MYNN/wrfout/'
time = pd.date_range('2017-'+str(current_month)+'-'+str(current_day)+' 9:00','2017-'+str(current_month)+'-'+str(current_day+1)+' 6:00',freq='15min')
wrf_dir = wrf_dir+time[0].strftime('%Y%m%d')+'/'
f = Dataset(wrf_dir+'wrfout_d01_2017'+time[0].strftime('%m%d')+'_09_00_00')
lat = f['XLAT'][0,:,:]
lon = f['XLONG'][0,:,:]
NKX_loc = wrf_m.find_nearest(lat,lon,32.85,-117.12)
time_pd = pd.DataFrame(np.arange(len(time)),index=time,columns=['index'])
t_index = time_pd.between_time('11:15','12:00').values[:,0]
z = getvar(f,'z',timeidx=t_index[0])[:,NKX_loc[0][0],NKX_loc[1][0]]

'''ED - heat flux'''
thetaL = wrf_m.get_thetaL_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)
ED_wpthetap, k3 = wrf_m.get_eddy_diffusivity_flux(f,z,t_index,NKX_loc[0][0],NKX_loc[1][0],thetaL,True)
ED_wpthetap, k3 = rho*cpAir*ED_wpthetap, k3
ql = wrf_m.get_qL_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)
ax4.plot(ED_wpthetap,z[0:-1],'--',color='saddlebrown')
'''ED - moisture flux'''
qt = wrf_m.get_qt_at_time_index(f,t_index,NKX_loc[0][0],NKX_loc[1][0],False)
ED_wpqtp, k4 = wrf_m.get_eddy_diffusivity_flux(f,z,t_index,NKX_loc[0][0],NKX_loc[1][0],qt,True)
ED_wpqtp, k4 = rho*Lv*ED_wpqtp, k4
ax5.plot(ED_wpqtp,z[0:-1],'--',color='saddlebrown')
#plt.savefig('/mnt/lab_45d1/database/Sc_group/WRF_39_RAP/EDMF_contribution/EDMF.png',dpi=200,bbox_inches='tight')