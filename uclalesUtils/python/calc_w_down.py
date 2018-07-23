#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 2018
Calculate coefficients for downdraft model
qt_downdraft(zi) = qt_avg(zi) + alpha*PBL_Flux/sigma_w(zi)
where alpha is a constant to be determined, sigma_w is empirically derived std of vertical velocity 
@author: elynn
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

def holtslag_Moeng_91(wstar,ustar,zstar,z):
    sigma_w = np.zeros_like(z)
    for i in range(len(sigma_w)):
        sigma_w[i] = 1.3*wstar*((ustar/wstar)**3+0.6*z[i]/zstar)**(1./3.)*(1.-z[i]/zstar)**(1./2.)
    return sigma_w

def stohl_etal_05(wstar,ustar,zstar,z):
    sigma_w = np.zeros_like(z)
    for i in range(len(sigma_w)):
        sigma_w[i] = (1.2*(1.-0.9*z[i]/zstar)*(z[i]/zstar)**(2./3.)+(1.8-1.4*(z[i]/zstar))*ustar**2)**(1./2.)
    return sigma_w
rf01_ps = Dataset('/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/rf01.ps.nc')
rf01_ts = Dataset('/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/rf01.ts.nc')
z = rf01_ps['zm'][:]
t_ps = rf01_ps['time'][:]
t_ts = rf01_ts['time'][:]
ind = np.where(t_ps[0] == t_ts)[0][0] #match first time index in ps to ts 
vertical_fluxes = pd.read_csv('/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/Vertical_fluxes_all_time.csv',index_col=0)
vertical_fluxes = vertical_fluxes[vertical_fluxes.columns[2::4]].mean(axis=1)
PBL_flux = vertical_fluxes.loc[0.9:1.05].min()#vertical_fluxes.loc[0.96:1].mean()
#1. Check how sigma_w scale with downdraft velocity
plume = pd.read_csv('/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/plume_output/zi2_bar_5percent_plume/DYCOMS_RF01_orig_vertical_profile_all_variables_hourly_avg.csv',index_col=0)
downdraft_w = plume['Downdraft w']
#need w* (Ghonima et al., 2017): 2.5*g/theta_v0 * integral{zb}{zi} <w'theta'>dz
wp_thetap = rf01_ps['tot_tw'][0,:]
zb = rf01_ts['zb'][ind]
zi = rf01_ts['zi2_bar'][ind]
wp_thetap_inCloud = np.interp(np.linspace(zb,zi,50),z,wp_thetap)
w_star = (2.5*9.8/289*np.trapz(wp_thetap_inCloud,np.linspace(zb,zi,50)))**(1./3.)
u_star = (rf01_ps['tot_uw'][0,0]**2 + rf01_ps['tot_vw'][0,0]**2) ** (1./4.)
sigma_w = w_star*1.3*((u_star/w_star)**3+0.6*z/zi)**(1./3.)*(1-z/zi)**(1./2.)
sigma_w_interp = np.interp(downdraft_w.index,z/zi,sigma_w)
mu_w = np.linspace(0,1,101)
diff = np.zeros_like(mu_w)
for i in range(len(diff)):
    diff[i] = np.nansum(np.abs(downdraft_w.as_matrix()/w_star+mu_w[i]*sigma_w_interp[::-1]/w_star)**2)
mu_w = mu_w[np.argmin(diff)]
plt.plot(downdraft_w.as_matrix()/w_star,downdraft_w.index,'-o',color='k',label='Downdraft w')
plt.plot(-sigma_w_interp*mu_w/w_star,1.-downdraft_w.index,color='b',label=r'$\sigma_w$ Holtslag and Moeng (1991)')

sigma_w_stohl = stohl_etal_05(w_star,u_star,zi,z)
sigma_w_stohl_interp = np.interp(downdraft_w.index,z/zi,sigma_w_stohl)
mu_w_stohl = np.linspace(2,3,101)
diff = np.zeros_like(mu_w_stohl)
for i in range(len(diff)):
    diff[i] = np.nansum(np.abs(downdraft_w.as_matrix()/w_star+mu_w_stohl[i]*sigma_w_stohl_interp/w_star)**2)
mu_w_stohl = mu_w_stohl[np.argmin(diff)]
plt.plot(-sigma_w_stohl_interp*mu_w_stohl/w_star,1.-downdraft_w.index,color='g',label=r'$\sigma_w$ Stohl et al. (2005)')
plt.legend(loc=0)
plt.xlabel(r'$w_{dd}/w_*, \sigma_w/w_*$')
plt.ylabel('z/z_*')
#plt.savefig('/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/downdraft_model/plume_velocity_with_sigma_w.png',dpi=200,bbox_inches='tight')
alpha = sigma_w_interp[0]/PBL_flux*(-9.875E-5)
#print 'theta_dd - theta_env= ', 