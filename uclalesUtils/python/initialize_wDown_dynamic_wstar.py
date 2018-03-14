#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 2018
Generate vertical w_down
@author: elynn
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import Sc_methods as Sc_m
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)
from oct2py import octave
octave.addpath('/home/elynn/Documents/github/Sc-utils/Soundings')
'''Commonly used constants'''
R = 287. #K/kgK
g = 9.81 #gravity [m/s^2]
cpAir = 1005. #specific heat of dry air [J/kgK]
rho_air = 1.2 #density of dry air [kg/m^3]
Lz = 1000. #depth of PBL height
Lv = 2.5e6
theta_ref = 290. #thetaV reference, set to 290K
'''End constants'''

'''CONFIGURE HERE'''
file_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02_no_mean_uv_uvmean0/'
case_prefix = 'rf01_nouv'
reducets_name = case_prefix+'.ts.nc'
reduceps_name = case_prefix+'.ps.nc'
LES = Dataset(file_dir+case_prefix+'_stitched.nc')
ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ps_ncfile = Dataset(file_dir+reduceps_name)
ts_time = ts_ncfile['time'][:] #time in the ts file
z = LES['zm'][:]
run_time_res = 60. #LES time resolution
run_startTime = 10860. #LES start time
run_endTime = 14400. #LES end time
z_scale = 'zi2_bar'
vertical_fluxes = pd.read_csv(file_dir+'plume_output/zi2_bar_5percent_plume/Vertical_fluxes_all_time.csv',index_col=0)
plume = pd.read_csv(file_dir+'plume_output/zi2_bar_5percent_plume/Vertical_profile_all_variables_hourly_avg.csv',index_col=0)
'''ENG CONFIGURATION'''

rflx = LES['rflx']
t = LES['t']
q = LES['q']
p = LES['p']
ql = LES['l']
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
deltaR = np.zeros_like(ts_LES)
delta_thetaL = np.zeros_like(ts_LES)
delta_physicalT = np.zeros_like(ts_LES)
delta_qt = np.zeros_like(ts_LES)
thetaL_cloud = np.zeros_like(ts_LES)
TL_cloud = np.zeros_like(ts_LES)
B0 = np.zeros_like(ts_LES)
Qinv = np.zeros_like(ts_LES)
Qcbl = np.zeros_like(ts_LES)
Tc = 289. #reference cloud temperature
delta_t = ts_LES[1]-ts_LES[0]
for t_index in range(1,len(ts_LES)):
    previous_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index-1])[0][0]]
    current_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index])[0][0]] #match time stamp in ts file
    deltaR[t_index] = Sc_m.get_deltaR(rflx[t_index,:,:,:])
    current_t = np.average(np.average(t[t_index,:,:,:],axis=0),axis=0) 
    delta_thetaL[t_index] = Sc_m.find_thetaL_jump(current_t,z/current_zstar)
    current_qt = np.average(np.average(q[t_index,:,:,:],axis=0),axis=0) 
    delta_qt[t_index] = Sc_m.find_qt_jump(current_qt,z/current_zstar)
    BL_idx = np.where((z/current_zstar)<=0.95)
    thetaL_cloud[t_index] = np.average(current_t[BL_idx])
    current_p = np.average(np.average(p[t_index,:,:,:],axis=0),axis=0)  
    current_ql = np.average(np.average(ql[t_index,:,:,:],axis=0),axis=0)
    current_physical_T = np.zeros_like(current_t)    
    for i in range(len(current_physical_T)):
        theta = current_t[i] + Lv/cpAir*current_ql[i]
        current_physical_T[i] = theta/(100000./current_p[i])**(R/cpAir)
    delta_physicalT[t_index] = octave.TMP_Inversion_Strength_Cal(current_physical_T,z/1000.,z[0])
    theta = current_t[BL_idx[0][-1]] + Lv/cpAir*current_ql[BL_idx[0][-1]]
    T = theta/(100000./current_p[BL_idx[0][-1]])**(R/cpAir)
    TL_cloud[t_index] = (T*cpAir + g*z[BL_idx[0][-1]] - Lv*current_ql[BL_idx[0][-1]])/cpAir
    B0[t_index] = deltaR[t_index]*g/(rho_air*cpAir*TL_cloud[t_index])
    inv_seg = np.where(np.diff(current_physical_T)>0.)[0]
    b_above = g * (current_physical_T[inv_seg[-1]] - Tc) / Tc
    b_below = g * (current_physical_T[inv_seg[0]] - Tc) / Tc
    Qinv[t_index] = np.abs((current_zstar - previous_zstar)/delta_t * (b_above - b_below))
    Qcbl[t_index] = -B0[t_index]-Qinv[t_index]
w_star = (-np.trapz(Qcbl[1:],ts_LES[1:])/(ts_LES[-1]-ts_LES[1])*current_zstar)**(1./3.)
deltaR = deltaR[1:].mean()
delta_T = delta_physicalT[1:].mean()
delta_qt = delta_qt[1:].mean()
TL_cloud = TL_cloud[1:].mean()
lamda = 15.
w_e = deltaR/rho_air/cpAir/delta_T * (0.175 + 0.39*w_star**2*Tc/lamda/g/delta_T)
print 'w*:', w_star, 'deltaR:', deltaR, 'deltaT:', delta_T, \
      'delta qt:', delta_qt, 'T cloud:', TL_cloud, 'w_e:', w_e
t_ps = ps_ncfile['time'][:]
t_ts = ts_ncfile['time'][:]
ind = np.where(t_ps[0] == t_ts)[0][0] #match first time index in ps to ts 
vertical_fluxes = vertical_fluxes[vertical_fluxes.columns[2::4]].mean(axis=1)
PBL_flux = vertical_fluxes.loc[0.9:1.05].min()
downdraft_w = plume['Downdraft w']
wp_thetap = ps_ncfile['tot_tw'][0,:]
zb = ts_ncfile['zb'][ind]
zi = ts_ncfile['zi2_bar'][ind]
wp_thetap_inCloud = np.interp(np.linspace(zb,zi,50),z,wp_thetap)
u_star = Sc_m.get_ustar(ps_ncfile,z)
BL_idx = np.argmin(np.abs(z/zi-1))
u_star = u_star[BL_idx]
print 'u*:', u_star
sigma_w = Sc_m.holtslag_Moeng_91(w_star,u_star,zi,z)
sigma_w_interp = np.interp(downdraft_w.index,z/zi,sigma_w)
#mu_w = np.linspace(0,3,101)
#diff = np.zeros_like(mu_w)
#for i in range(len(diff)):
#    diff[i] = np.nansum(np.abs(downdraft_w.as_matrix()/w_star+mu_w[i]*sigma_w_interp[::-1]/w_star)**2)
#mu_w = mu_w[np.argmin(diff)]
mu_w = 1.8 #derived from DYCOMS RF01 regular
plt.plot(downdraft_w.as_matrix()/w_star,downdraft_w.index,'-o',color='k',label='Downdraft w')
plt.plot(-sigma_w_interp*mu_w/w_star,1.-downdraft_w.index,color='b',label=r'$\sigma_w$ Holtslag and Moeng (1991)')
plt.legend(loc=0)
plt.xlabel(r'$w_{dd}/w_*, \sigma_w/w_*$')
plt.ylabel(r'$z/z_*$')
plt.xlim([-1.3,0])