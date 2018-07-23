#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:07:42 2018

@author: elynn
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.1)
from oct2py import octave
octave.addpath('/home/elynn/Documents/github/Sc-utils/Soundings')
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
R = 287. #K/kgK
g = 9.81 #gravity [m/s^2]
cpAir = 1005. #specific heat of dry air [J/kgK]
rho_air = 1.2 #density of dry air [kg/m^3]
Lz = 1000. #depth of PBL height
Lv = 2.5e6
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
p = LES['p']
ql = LES['l']
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
deltaR = np.zeros_like(ts_LES)
delta_thetaL = np.zeros_like(ts_LES)
delta_physicalT = np.zeros_like(ts_LES)
delta_qt = np.zeros_like(ts_LES)
thetaL_cloud = np.zeros_like(ts_LES)
TL_cloud = np.zeros_like(ts_LES)
for t_index in range(len(ts_LES)):
    current_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index])[0][0]] #match time stamp in ts file
    current_R = np.average(np.average(rflx[t_index,:,:,:],axis=0),axis=0)
    deltaR[t_index] = current_R.max() - current_R.min()
    current_t = np.average(np.average(t[t_index,:,:,:],axis=0),axis=0) 
    delta_thetaL[t_index] = find_thetaL_jump(current_t,z/current_zstar)
    current_qt = np.average(np.average(q[t_index,:,:,:],axis=0),axis=0) 
    delta_qt[t_index] = find_qt_jump(current_qt,z/current_zstar)
    BL_idx = np.where((z/current_zstar)<=0.95)
    thetaL_cloud[t_index] = np.average(current_t[BL_idx])
    current_p = np.average(np.average(p[t_index,:,:,:],axis=0),axis=0)  
    current_ql = np.average(np.average(ql[t_index,:,:,:],axis=0),axis=0)
    ###Call octave to compute delta T jump using Xiaohui's algorithm
    current_physical_T = np.zeros_like(current_t)    
    for i in range(len(current_physical_T)):
        theta = current_t[i] + Lv/cpAir*current_ql[i]
        current_physical_T[i] = theta/(100000./current_p[i])**(R/cpAir)
    delta_physicalT[t_index] = octave.TMP_Inversion_Strength_Cal(current_physical_T,z/1000.,z[0])
    theta = current_t[BL_idx[0][-1]] + Lv/cpAir*current_ql[BL_idx[0][-1]]
    T = theta/(100000./current_p[BL_idx[0][-1]])**(R/cpAir)
    TL_cloud[t_index] = (T*cpAir + g*z[BL_idx[0][-1]] - Lv*current_ql[BL_idx[0][-1]])/cpAir
deltaR = deltaR.mean()
delta_T = delta_physicalT.mean()
delta_qt = delta_qt.mean()
TL_cloud = TL_cloud.mean()
w_star = 1.0 #approximation for now, needs to be
lamda = 15.
w_e = deltaR/rho_air/cpAir/delta_T * (0.175 + 0.39*w_star**2*TL_cloud/lamda/g/delta_T)

thetaL_diff = -0.022 #calculated from regular DYCOMS, 5% plume, z/zi between [0,0.95]
qt_diff = -9.875E-5 #calculated from regular DYCOMS, 5% plume, z/zi between [0,0.95]

alpha1 = qt_diff*w_star/w_e/delta_qt
alpha2 = -thetaL_diff*w_star/deltaR/g*rho_air*cpAir*TL_cloud