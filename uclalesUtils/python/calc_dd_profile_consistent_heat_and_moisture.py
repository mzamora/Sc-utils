#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 2018
Calculate coefficients for downdraft model (heat + moisture)
Following Kohler's ECMWF notes
This is to test if the coefficient for heat and moisture is consistent
theta_diff = b*(deltaR+wPrime_thetaPrime_ent)/w*
qt_diff = b*(wPrime_qtPrime)/w*
w* = (g/thetaV*(deltaR+wPrime_thetaPrime_ent)*Lz)^(1/3)
@author: elynn
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.1)
'''CONFIGURE HERE'''
file_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/'
reducets_name = 'rf01.ts.nc'
LES = Dataset(file_dir+'rf01_stitched.nc')
ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
z = LES['zm'][:]
run_time_res = 60. #LES time resolution
run_startTime = 10860. #LES start time
run_endTime = 14400. #LES end time
z_scale = 'zi2_bar'
'''ENG CONFIGURATION'''
'''Commonly used constants'''
g = 9.81 #gravity [m/s^2]
cpAir = 1005. #specific heat of dry air [J/kgK]
rho_air = 1.2 #density of dry air [kg/m^3]
Lz = 1000. #depth of PBL height
theta_ref = 290. #thetaV reference, set to 290K
'''End constants'''

def find_thetaL_jump(thetaL,z_normalized):
    '''Find thetaL inversion jump
    Parameters
    ----------
    thetaL: np.array
        Liquid water potential temperature
    z_normalized: np.array
        Height normalized by inversion base
    Returns
    -------
    delta_thetaL: float
        Inversion jump [K]
    '''
    slope = np.diff(thetaL)
    max_ind = np.argmax(slope) #find the maximum index first, then starts iterating
    while slope[max_ind+1]<slope[max_ind]: #starts iterating up until local minimum
        max_ind = max_ind + 1
    BL_ind = np.where(z_normalized<=1)[0]
    BL_thetaL = np.average(thetaL[BL_ind])
    delta_thetaL = thetaL[max_ind] - BL_thetaL
    return delta_thetaL

def find_qt_jump(qt,z_normalized):
    '''Find qt inversion jump
    Parameters
    ----------
    qt: np.array
        Total water mixing ratio
    z_normalized: np.array
        Height normalized by inversion base
    Returns
    -------
    delta_qt: float
        Inversion jump [K]
    '''
    slope = np.diff(qt)
    max_ind = np.argmin(slope) #find the maximum index first, then starts iterating
    while slope[max_ind+1]>slope[max_ind]: #starts iterating up until local minimum
        max_ind = max_ind + 1
    BL_ind = np.where(z_normalized<=1)[0]
    BL_qt = np.average(qt[BL_ind])
    delta_qt = qt[max_ind] - BL_qt
    return delta_qt

#1. Get buoyancy flux - radiative cooling term from LES (rflx diff.)
rflx = LES['rflx']
t = LES['t']
q = LES['q']
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
deltaR = np.zeros_like(ts_LES)
delta_thetaL = np.zeros_like(ts_LES)
delta_qt = np.zeros_like(ts_LES)
for t_index in range(len(ts_LES)):
    current_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index])[0][0]] #match time stamp in ts file
    current_R = np.average(np.average(rflx[t_index,:,:,:],axis=0),axis=0)
    deltaR[t_index] = current_R.max() - current_R.min()
    current_t = np.average(np.average(t[t_index,:,:,:],axis=0),axis=0) 
    delta_thetaL[t_index] = find_thetaL_jump(current_t,z/current_zstar)
    current_qt = np.average(np.average(q[t_index,:,:,:],axis=0),axis=0) 
    delta_qt[t_index] = find_qt_jump(current_qt,z/current_zstar)
deltaR = deltaR.mean()
delta_thetaL = 8.#delta_thetaL.mean() #8.
delta_qt = delta_qt.mean()

#2. Get buoyancy flux - entrainment part
####IMPORTANT!!! For the time being, assume entrainment part is 30% of radiaitve cooling
####This assumption is consistent with Kohler's value as well as Davini et al., 2017 paper
####Needs to be updated when officially implemented (e.g. get w_e from Ghonima et al., 2017 eq.3)
wPrime_thetaPrime_ent = 0.3*deltaR
w_e = wPrime_thetaPrime_ent/cpAir/rho_air/delta_thetaL

#3. Get w_star
w_star = (g/theta_ref*(deltaR-wPrime_thetaPrime_ent)/cpAir/rho_air*Lz)**(1./3.)

#4. Get coefficient b for heat 
thetaL_diff = -0.022 #calculated from regular DYCOMS, 5% plume, z/zi between [0,0.95]
b_heat = -thetaL_diff/(deltaR-wPrime_thetaPrime_ent)*cpAir*rho_air*w_star

#5. Get coefficient b for moisture
qt_diff = -9.875E-5 #calculated from regular DYCOMS, 5% plume, z/zi between [0,0.95]
b_moisture = qt_diff*w_star/w_e/delta_qt

print 'b_heat:', b_heat, 'b_moisture:', b_moisture