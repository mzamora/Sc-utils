#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 2018
@author: elynn

Get Richardson number at inversion
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

'''Configure here'''
file_dir = '/home/elynn/Documents/uclales/drycbl/hr_3_4_1min/'
reducets_name = 'dcbl_x.ts.nc'
run_stitched_file = 'dcbl_x_hr34_1min_stitched.nc'
run_time_res = 60. #LES time resolution [s]
run_startTime = 10860. #LES start time
run_endTime = 14400. #LES end time
z_scale = 'zi2_bar'
output_dir = '/home/elynn/Documents/uclales/drycbl/hr_3_4_1min/'
'''End configuration'''

'''Read time statistic file and stitched LES file'''
ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
les_ncfile = Dataset(file_dir+run_stitched_file) #stitched LES output
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
u = les_ncfile['u']
v = les_ncfile['v']
z = les_ncfile['zm'][:]
q = les_ncfile['q']
#ql = les_ncfile['l']
theta = les_ncfile['t']
Lv = 2.5E6
cpAir = 1005.
#theta = thetaL[:,:,:,:]+Lv/cpAir*ql[:,:,:,:]
thetaV = theta*(1.0+0.61*q[:,:,:,:])
Ri = np.zeros(len(ts_LES))
for t_index in range(len(ts_LES)):
    current_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index])[0][0]] #match time stamp in ts file
    zstar = z/current_zstar
    ind = np.argmin(np.abs(zstar-1)) #find the closest index to PBL height
    current_thetaVL = np.average(np.average(thetaV[t_index,:,:,:],axis=0),axis=0)
    current_u = np.average(np.average(u[t_index,:,:,:],axis=0),axis=0)
    current_v = np.average(np.average(v[t_index,:,:,:],axis=0),axis=0)
    du = current_u[ind+1]-current_u[ind-1]
    dv = current_v[ind+1]-current_v[ind-1]
    Ri[t_index] = 9.8*(current_thetaVL[ind+1]-current_thetaVL[ind-1])*np.diff(z)[ind]/current_thetaVL[ind]/(du**2+dv**2)
            
np.savetxt(output_dir+'Richardson_number_ts.csv',Ri,delimiter=',')