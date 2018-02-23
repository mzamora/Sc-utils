#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 2018
@author: elynn

Get vertical fluxes
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

'''Configure here'''
file_dir = '/home/elynn/Documents/uclales/cgils_s12_ctl/hr_3_4_1min/'
reducets_name = 'cgils_s12_ctl.ts.nc'
run_stitched_file = 'cgils_s12_ctl_hr_3_4_1min_stitched.nc'
run_time_res = 60. #LES time resolution
run_startTime = 10860. #LES start time
run_endTime = 14400. #LES end time
z_scale = 'zi2_bar'
output_dir = '/home/elynn/Documents/uclales/cgils_s12_ctl/hr_3_4_1min/plume_output/zi2_bar_5percent_plume/'
'''End configuration'''

ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
les_ncfile = Dataset(file_dir+run_stitched_file) #stitched LES output
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
x = les_ncfile['xm'][:]
y = les_ncfile['ym'][:]
z = les_ncfile['zm'][:]
w_bar = np.average(les_ncfile['w'][:,:,:,:],axis=0)
thetaL_bar = np.average(les_ncfile['t'][:,:,:,:],axis=0)
thetaV = les_ncfile['t'][:,:,:,:]*(1.0+0.608*les_ncfile['q'][:,:,:,:])
thetaV_bar = np.average(thetaV,axis=0)
qT_bar = np.average(les_ncfile['q'][:,:,:,:],axis=0)
ql_bar = np.average(les_ncfile['l'][:,:,:,:],axis=0)
flag = True
for t_index in range(len(ts_LES)):
    wp_thetaVp = np.zeros(len(z)) #w'tv'
    wp_thetaLp = np.zeros(len(z)) #w'tl'
    wp_qtp = np.zeros(len(z)) #w'qt'
    wp_qlp = np.zeros(len(z)) #w'ql'
    current_zstar = ts_ncfile[z_scale][np.where(ts_time == ts_LES[t_index])[0][0]] #match time stamp in ts file
    for k in range(len(z)):
        wp_thetaVp_k = np.average((les_ncfile['w'][t_index,:,:,k]-w_bar[:,:,k]) * (thetaV[t_index,:,:,k]-thetaV_bar[:,:,k]))
        wp_thetaVp[k] = wp_thetaVp_k
        wp_thetaLp_k = np.average((les_ncfile['w'][t_index,:,:,k]-w_bar[:,:,k]) * (les_ncfile['t'][t_index,:,:,k]-thetaL_bar[:,:,k]))
        wp_thetaLp[k] = wp_thetaLp_k
        wp_qtp_k = np.average((les_ncfile['w'][t_index,:,:,k]-w_bar[:,:,k]) * (les_ncfile['q'][t_index,:,:,k]-qT_bar[:,:,k]))
        wp_qtp[k] = wp_qtp_k        
        wp_qlp_k = np.average((les_ncfile['w'][t_index,:,:,k]-w_bar[:,:,k]) * (les_ncfile['l'][t_index,:,:,k]-ql_bar[:,:,k]))
        wp_qlp[k] = wp_qlp_k
    zcoord = z/current_zstar #current z scaled by current PBL height
#    wp_thetaVp = wp_thetaVp/wp_thetaVp[0] #scaled by sfc value
#    wp_thetaLp = wp_thetaLp/wp_thetaLp[0] #scaled by sfc value
#    wp_qtp = wp_qtp/wp_qtp[0] #scaled by sfc value
#    wp_qlp = wp_qlp/wp_qlp[0] #scaled by sfc value
    wp_thetaVp = np.interp(np.arange(0,1.21,0.01),zcoord,wp_thetaVp)
    wp_thetaLp = np.interp(np.arange(0,1.21,0.01),zcoord,wp_thetaLp)
    wp_qtp = np.interp(np.arange(0,1.21,0.01),zcoord,wp_qtp)
    wp_qlp = np.interp(np.arange(0,1.21,0.01),zcoord,wp_qlp)
    if flag:
        output = pd.DataFrame(wp_thetaVp,index=np.arange(0,1.21,0.01),columns=['wp_thetaVp'])
        output['wp_thetaLp'] = wp_thetaLp
        output['wp_qtp'] = wp_qtp
        output['wp_qlp'] = wp_qlp
        flag = False
    else:
        current = pd.DataFrame(wp_thetaVp,index=np.arange(0,1.21,0.01),columns=['wp_thetaVp'])
        current['wp_thetaLp'] = wp_thetaLp
        current['wp_qtp'] = wp_qtp
        current['wp_qlp'] = wp_qlp
        output = pd.concat((output,current),axis=1)
output.to_csv(output_dir+'Vertical_fluxes_all_time.csv')
