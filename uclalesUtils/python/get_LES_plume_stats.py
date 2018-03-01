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
z_scale = 'zi2_bar'
output_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/plume_output/zi2_bar_5percent_plume/'
'''End configuration'''

'''Read time statistic file and stitched LES file'''
ts_ncfile = Dataset(file_dir+reducets_name) #time statistics from reducets
ts_time = ts_ncfile['time'][:] #time in the ts file
les_ncfile = Dataset(file_dir+run_stitched_file) #stitched LES output
ts_LES = np.arange(run_startTime,run_endTime+run_time_res,run_time_res) #input LES time stamp
plume_vars = ['t','q','l','w'] #plume profiles: liquid potential temperature, total water mixing ratio, liquid water mixing ratio, vertical velocity
#plume_vars = ['t','q','w']
w = les_ncfile['w']
t = les_ncfile['t']
z = les_ncfile['zm'][:]

'''3 loops - time, variables, vertical slice'''
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
                selector = np.where(current_w<=np.percentile(current_w,7)) #find downdraft
                current = current_var[selector].mean()
                vertical_var.append(current)
                selector = np.where(current_w>=np.percentile(current_w,93)) #find updraft
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
        domain_var = np.interp(np.arange(0,1.21,0.01),zcoord,domain_var)
        vertical_var = np.interp(np.arange(0,1.21,0.01),zcoord,vertical_var)
        vertical_var_ud = np.interp(np.arange(0,1.21,0.01),zcoord,vertical_var_ud)
        if flag:
            output = pd.DataFrame(domain_var,index=np.arange(0,1.21,0.01),columns=['domain_'+plume_vars[vc]])
            output['updraft_'+plume_vars[vc]] = vertical_var_ud
            output['downdraft_'+plume_vars[vc]] = vertical_var
            flag = False
        else:
            output['domain_'+plume_vars[vc]] = domain_var
            output['updraft_'+plume_vars[vc]] = vertical_var_ud
            output['downdraft_'+plume_vars[vc]] = vertical_var
    output.to_csv(output_dir+'Vertical_profile_at_time_'+str(t_index+1)+'.csv')

'''plotting function'''
def plot_hourly_avg_vert_profile(plume_vars,output_dir):
    xlabels = [r'$\theta_l$ [K]',r'$q_T$ [g/kg]',r'$q_l$ [g/kg]',r'$w$ [m/s]']
#    xlabels = [r'$\theta_l$ [K]',r'$q_T$ [g/kg]',r'$w$ [m/s]']
    xlimits = [[289.1,289.5],[8.5,9.5],[-0.01,0.7],[-2.0,2.0]] #DYCOMS RF01
#    xlimits = [[288.2,288.7],[9.0,10.0],[-0.01,0.7],[-2.0,2.0]] #CGILS S12 ctl
#    xlimits = [[303.2,304.2],[9.8,10.2],[-2.2,2.2]] #DryCBL
    flag = True
    for t_index in range(1,61):
        current = pd.read_csv(output_dir+'Vertical_profile_at_time_'+str(t_index)+'.csv',index_col=0) 
        if flag:
            output = current
            flag = False
        else:
            output = pd.concat([output,current],axis=1)
    plt.figure(figsize=(15,5))
    for i in range(len(plume_vars)):
        plume_var = plume_vars[i]
        if (plume_var == 'q') or (plume_var == 'l'):
            domain = output['domain_'+plume_var].mean(axis=1)*1000.
            updraft = output['updraft_'+plume_var].mean(axis=1)*1000.
            downdraft = output['downdraft_'+plume_var].mean(axis=1)*1000.  
        else:
            domain = output['domain_'+plume_var].mean(axis=1)
            updraft = output['updraft_'+plume_var].mean(axis=1)
            downdraft = output['downdraft_'+plume_var].mean(axis=1)     
        plt.subplot(1,4,i+1)
        plt.plot(domain.iloc[:],domain.index[:],color='k',label='domain')
        plt.plot(updraft.iloc[:],domain.index[:],color='r',label='updraft')
        plt.plot(downdraft.iloc[:],domain.index[:],color='b',label='downdraft')
        plt.xlabel(xlabels[i])
        plt.xlim(xlimits[i])
        plt.ylabel(r'$z/z_*$ [-]')
        if i==0:
            plt.legend(loc=0) 
    suptitle = plt.suptitle('Hour3-4 average', y=1.05)
    plt.tight_layout()
#    plt.savefig(output_dir+'Plume_vars_hour_34_avg.png',dpi=200,bbox_inches='tight',bbox_extra_artists=[suptitle])
    return output
output = plot_hourly_avg_vert_profile(plume_vars,output_dir)
#dd_t = output['downdraft_t'].mean(axis=1)
#domain_t = output['domain_t'].mean(axis=1)
#vert = pd.concat((dd_t,domain_t),axis=1)
#vert = pd.concat((dd_t,domain_t),axis=1,names=['Downdraft','Domain'])
#vert.to_csv(output_dir+'Vertical_profile_hourly_avg.csv')
dd_q = output['downdraft_q'].mean(axis=1)
ud_q = output['updraft_q'].mean(axis=1)
domain_q = output['domain_q'].mean(axis=1)
dd_l = output['downdraft_l'].mean(axis=1)
ud_l = output['updraft_l'].mean(axis=1)
domain_l = output['domain_l'].mean(axis=1)
dd_t = output['downdraft_t'].mean(axis=1)
ud_t = output['updraft_t'].mean(axis=1)
domain_t = output['domain_t'].mean(axis=1)
dd_w = output['downdraft_w'].mean(axis=1)
ud_w = output['updraft_w'].mean(axis=1)
vert = pd.concat((dd_t,ud_t,domain_t,dd_q,ud_q,domain_q,dd_l,ud_l,domain_l,dd_w,ud_w),axis=1)
vert.columns = ['Downdraft thetaL','Updraft thetaL','Domain thetaL',\
                'Downdraft qt','Updraft qt','Domain qt',\
                'Downdraft ql','Updraft ql','Domain ql',\
                'Downdraft w','Updraft w']
vert.to_csv(output_dir+'Vertical_profile_all_variables_hourly_avg.csv')