#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 2018
@author: elynn

Get entrainment velocity from UCLALES output
Mass-based entrainment velocity [m/s]:
    w_e = dzi_dt + vDotGradzi - w(IBH)
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

'''Configure here'''
file_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/'
reducets_name = 'rf01.ts.nc'
run_stitched_file = 'rf01_stitched.nc'
run_time_res = 60. #LES time resolution
run_startTime = 10860. #LES start time
run_endTime = 14400. #LES end time
z_scale = 'zi2_bar'
output_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/'
'''End configuration'''

ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
les_ncfile = Dataset(file_dir+run_stitched_file) #stitched LES output
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
x = les_ncfile['xm'][:]
y = les_ncfile['ym'][:]
z = les_ncfile['zm'][:]
ql = les_ncfile['l']
dt = run_time_res
h = np.diff(x)[0]
PBL_ind = np.zeros((len(ts_LES),len(y),len(x)),dtype=int)
for t_index in range(len(ts_LES)):
    for i in range(len(y)):
        for j in range(len(x)):
            current_ql = ql[t_index,i,j,:]
            if np.max(current_ql)>0: #a cloud exist
                cloud_top_index = np.nonzero(current_ql)[0][-1]
                PBL_ind[t_index,i,j] = cloud_top_index
            else:
                PBL_ind[t_index,i,j] = -999            
    median = np.median(PBL_ind[t_index,:,:])
    PBL_ind[PBL_ind<0] = int(median) #grid points without clouds will be filled with median PBL height
div = 3.75e-6    
w_e = np.zeros((len(ts_LES),len(y),len(x)))
for t_index in range(len(ts_LES)):
    for i in range(1,len(y)-1):
        for j in range(1,len(x)-1):
            if t_index == 0: #forward difference at the start edge
                dzi_dt = (1./dt) * (z[PBL_ind[t_index+1,i,j]] - z[PBL_ind[t_index,i,j]])
            elif t_index == len(ts_LES)-1: #backward difference at the end edge
                dzi_dt = (1./dt) * (z[PBL_ind[t_index,i,j]] - z[PBL_ind[t_index-1,i,j]])
            else: #otherwise, 2nd order central difference for points in the middle
                dzi_dt = (1./(2.*dt)) * (z[PBL_ind[t_index+1,i,j]] - z[PBL_ind[t_index-1,i,j]])
            
            meanU = np.average(les_ncfile['u'][t_index,j,i,0:PBL_ind[t_index,i,j]+1])
            meanV = np.average(les_ncfile['v'][t_index,j,i,0:PBL_ind[t_index,i,j]+1])
            dzi_dx = (1./(2.*h)) * (z[PBL_ind[t_index,i,j+1]] - z[PBL_ind[t_index,i,j-1]])
            dzi_dy = (1./(2.*h)) * (z[PBL_ind[t_index,i+1,j]] - z[PBL_ind[t_index,i-1,j]])
            vDotGradzi = meanU*dzi_dx + meanV*dzi_dy
            w_e[t_index,i,j] = dzi_dt + vDotGradzi + div * z[PBL_ind[t_index,i,j]]          
w_e_avg = np.average(w_e,axis=1)
w_e_avg = np.average(w_e_avg,axis=1)
w_e_avg = pd.DataFrame(w_e_avg*1000.,index=ts_LES,columns=['Entrainment [mm/s]'])
w_e_avg.to_csv(output_dir+'Entrainment_rate.csv')