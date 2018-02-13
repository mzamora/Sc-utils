#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 2018
@author: elynn

Get updraft/downdraft statistics from UCLALES output
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
z_scale = 'zi3_bar'
output_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/plume_output/zi3_bar_5percent_plume/'
'''End configuration'''

'''Read time statistic file and stitched LES file'''
ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
les_ncfile = Dataset(file_dir+run_stitched_file) #stitched LES output
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
plume_vars = ['t','q','l','w'] #plume profiles: liquid potential temperature, total water mixing ratio, liquid water mixing ratio, vertical velocity
w = les_ncfile['w']
ql = les_ncfile['l']
t = les_ncfile['t']
z = les_ncfile['zm'][:]
for t_index in range(60):
    flag = True #first time in each t_index
    for vc in range(len(plume_vars)): #loop over plume profiles
        current_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index])[0][0]] #match time stamp in ts file
        var = les_ncfile[plume_vars[vc]]
        vertical_var = []
        vertical_var_ud = []
        domain_var = []
        for i in range(len(z)): #loop over the vertical slice 
            current_var = var[t_index,:,:,i]
            current_w = w[t_index,:,:,i]
            if np.abs(current_w).max()>0:
                selector = np.where(current_w<=np.percentile(current_w,5)) #find downdraft
                current = current_var[selector].mean()
                vertical_var.append(current)
                selector = np.where(current_w>=np.percentile(current_w,95)) #find updraft
                current = current_var[selector].mean()
                vertical_var_ud.append(current)
            else:
                vertical_var.append(np.nan)
                vertical_var_ud.append(np.nan)
            domain_var.append(current_var.mean())
        vertical_var = np.array(vertical_var) #list to array
        vertical_var_ud = np.array(vertical_var_ud) #list to array
        domain_var = np.array(domain_var) #list to array
        zcoord = z/current_zstar #current z scaled by current PBL height
        domain_var = np.interp(np.arange(0,1.01,0.01),zcoord,domain_var)
        vertical_var = np.interp(np.arange(0,1.01,0.01),zcoord,vertical_var)
        vertical_var_ud = np.interp(np.arange(0,1.01,0.01),zcoord,vertical_var_ud)
        if flag:
            output = pd.DataFrame(domain_var,index=np.arange(0,1.01,0.01),columns=['domain_'+plume_vars[vc]])
            output['updraft_'+plume_vars[vc]] = vertical_var_ud
            output['downdraft_'+plume_vars[vc]] = vertical_var
            flag = False
        else:
            output['domain_'+plume_vars[vc]] = domain_var
            output['updraft_'+plume_vars[vc]] = vertical_var_ud
            output['downdraft_'+plume_vars[vc]] = vertical_var
    output.to_csv(output_dir+'Vertical_profile_at_time_'+str(t_index+1)+'.csv')