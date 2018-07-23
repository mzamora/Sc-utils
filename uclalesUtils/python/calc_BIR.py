#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 2018
Calculate buoyancy flux integral ratio (BIR) following: 
Stevens, B., 2000. Cloud transitions and decoupling in shear-free stratocumulus-topped boundary layers 27, 2557–2560.
@author: elynn
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import Sc_methods as Sc_m
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

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
file_dir = '/home/elynn/Documents/uclales/NKX_mz/hr_3_4_1min/'
case_prefix = 'nkxfullradlong2'
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
'''ENG CONFIGURATION'''
wp_thetaVp = vertical_fluxes[vertical_fluxes.columns[0::4]].mean(axis=1).loc[0:1]
z = np.array(wp_thetaVp.index.tolist())*ts_ncfile[z_scale][:].mean()
wp_thetaVp = wp_thetaVp.as_matrix()
negative = np.where(wp_thetaVp<0)
positive = np.where(wp_thetaVp>0)
num = np.trapz(wp_thetaVp[negative[0]],z[negative[0]])
denum = np.trapz(wp_thetaVp[positive[0]],z[positive[0]])
BIR = -num/denum
print BIR*100.
