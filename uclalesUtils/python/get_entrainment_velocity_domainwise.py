#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 5 2018
@author: elynn

Get entrainment velocity from UCLALES output
Mass-based entrainment velocity [m/s]:
    w_e = dzi_dt + vDotGradzi - Dw(IBH)
Domain-wise: periodic condition, skip over vDotGradzi
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

'''Configure here'''
file_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02_no_mean_uv_uvmean0/hr_3_4_1min/'
reducets_name = 'rf01_nouv.ts.nc'
run_stitched_file = 'rf01_nouv_stitched.nc'
run_time_res = 60. #LES time resolution
run_startTime = 10860. #LES start time
run_endTime = 14400. #LES end time
z_scale = 'zi2_bar'
output_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02_no_mean_uv_uvmean0/hr_3_4_1min/'
'''End configuration'''
ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
zi = ts_ncfile[z_scale][:]
dt = np.diff(ts_time)[0]
w_e = np.zeros_like(ts_time)
div = 3.75e-6  #5.18016507e-6#CGILS#
for t_index in range(len(ts_time)):
    if t_index == 0: #forward difference at the start edge
        dzi_dt = (1./dt) * (zi[t_index+1] - zi[t_index])
    elif t_index == len(ts_time)-1: #backward difference at the end edge
        dzi_dt = (1./dt) * (zi[t_index] - zi[t_index-1])
    else: #otherwise, 2nd order central difference for points in the middle
        dzi_dt = (1./(2.*dt)) * (zi[t_index+1] - zi[t_index-1])    
    w_e[t_index] = dzi_dt + div * zi[t_index]
print w_e.mean()