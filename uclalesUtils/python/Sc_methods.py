import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

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
    max_ind = np.argmax(slope) #find the m/home/elynn/Documents/uclales/dycomsrf01_br02_no_mean_uv_uvmean0/hr_3_4_1min/aximum index first, then starts iterating
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

def holtslag_Moeng_91(wstar,ustar,zstar,z):
	'''Sigma w used in Siebesma et al., 2007
	'''
	sigma_w = np.zeros_like(z)
	for i in range(len(sigma_w)):
	    sigma_w[i] = 1.3*wstar*((ustar/wstar)**3+0.6*z[i]/zstar)**(1.0/3.0)*(1.0-z[i]/zstar)**(1.0/2.0)
	return sigma_w

def get_deltaR(rflx,radtype,ps_ncfile,zi,zb,z,inv_top):
	'''
	Get total radiative flux jump
	Parameters
	----------
	rflx: np.array (3d), z in last index
		Total radiative flux in (x, y, z)
	radtype: int
		Radiation scheme in uclales, 4 for full radiation, 2 for DYCOMS (simplified)
	pc_ncfile: opened nc file
		Profile netcdf file
	inv_top: int
		Inversion top index
	Returns
	-------
	delta_R: float
		Total radiative flux jump
	'''
	if radtype == 2:
		current_R = np.average(np.average(rflx[:,:,:],axis=0),axis=0)
		delta_R = current_R.max() - current_R.min()
	elif radtype == 4:
		rad = ps_ncfile['rflx'][:] - ps_ncfile['sflx'][:]
		if rad.shape[0]>1:
			rad = np.average(rad[1:,:],axis=0)
		else:
			rad = rad[0,:]
		BL_idx = np.argmin(np.abs(z/zi-1))
		zb_idx = np.argmin(np.abs(z/zb-1))
		delta_R = rad[inv_top] - rad[zb_idx]#rad.min()
	return delta_R	

def get_ustar(ps_ncfile,z):
	'''
	Get u* as formulated (upwp^2+vpwp^2)
	Parameters
	----------
	ps_nc: ps netcdf file
		Opened nc file
	z: np.array
		Vertical grid
	Returns
	-------
	u_star: np.array like z
		u* in the vertical grid
	'''
	u_star = np.zeros_like(z)
	upwp = np.average(ps_ncfile['tot_uw'][:,:],axis=0)
	vpwp = np.average(ps_ncfile['tot_vw'][:,:],axis=0)
	for i in range(len(u_star)):
	    u_star[i] = (upwp[i]**2 + vpwp[i]**2) ** (1./4.)
	return u_star