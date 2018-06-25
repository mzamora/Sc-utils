"""

Methods for post-processing UCSD WRF data

@author: Elynn
"""
import numpy as np
import datetime
import pytz
import pandas as pd
import os
from wrf import getvar, destagger
import fnmatch
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from oct2py import octave
from math import radians, cos, sin, asin, sqrt
rho = 1.28
Lv = 2.5e6
cpAir = 1005. #specific heat of dry air [J/kgK]
def get_thetaL_at_time_index(f,t_index,loci,locj,average_flag):
    '''Get liquid potential temeprature at the specified time index
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    t_index: int, or array of ints
        index/indices for time
    loci, locj: int
        index to lat and lon
    average_flag: boolean
        whether to avearage time indicies when len(t_index)>1
    Returns
    -------
    thetaL: np.array of float
        Liquid water potential temperature
    '''
    if type(t_index) is int: #only 1 single time index is specified
        theta = f['T'][t_index,:,loci,locj]+300.
        ql = f['QCLOUD'][t_index,:,loci,locj]
        thetaL = theta - Lv/cpAir * ql
        average_flag = False #no average allowed because only 1 t_index specified
    else: #t_index is a list or array
        theta = np.zeros((len(t_index),f['T'].shape[1]))
        ql = np.zeros((len(t_index),f['T'].shape[1]))
        thetaL = np.zeros((len(t_index),f['T'].shape[1]))
        for i in range(len(t_index)):
            theta = f['T'][t_index[i],:,loci,locj]+300.
            ql = f['QCLOUD'][t_index[i],:,loci,locj]
            thetaL[i,:] = theta - Lv/cpAir * ql
    if average_flag:
        thetaL = np.average(thetaL,axis=0)
    return thetaL

def get_qL_at_time_index(f,t_index,loci,locj,average_flag):
    '''Get liquid water mixing ratio at the specified time index
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    t_index: int, or array of ints
        index/indices for time
    loci, locj: int
        index to lat and lon
    average_flag: boolean
        whether to avearage time indicies when len(t_index)>1
    Returns
    -------
    thetaL: np.array of float
        Liquid water potential temperature
    '''
    if type(t_index) is int: #only 1 single time index is specified
        ql = f['QCLOUD'][t_index,:,loci,locj]
        average_flag = False #no average allowed because only 1 t_index specified
    else: #t_index is a list or array
        ql = np.zeros((len(t_index),f['T'].shape[1]))
        thetaL = np.zeros((len(t_index),f['T'].shape[1]))
        for i in range(len(t_index)):
            ql[i,:] = f['QCLOUD'][t_index[i],:,loci,locj]
    if average_flag:
        thetaL = np.average(thetaL,axis=0)
    return ql

def get_qt_at_time_index(f,t_index,loci,locj,average_flag):
    '''Get total water mixing ratio at the specified time index
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    t_index: int, or array of ints
        index/indices for time
    loci, locj: int
        index to lat and lon
    average_flag: boolean
        whether to avearage time indicies when len(t_index)>1
    Returns
    -------
    thetaL: np.array of float
        Liquid water potential temperature
    '''
    if type(t_index) is int: #only 1 single time index is specified
        ql = f['QCLOUD'][t_index,:,loci,locj]
        qv = f['QVAPOR'][t_index,:,loci,locj]
        average_flag = False #no average allowed because only 1 t_index specified
    else: #t_index is a list or array
        ql = np.zeros((len(t_index),f['T'].shape[1]))
        qv = np.zeros((len(t_index),f['T'].shape[1]))
        for i in range(len(t_index)):
            ql[i,:] = f['QCLOUD'][t_index[i],:,loci,locj]
            qv[i,:] = f['QVAPOR'][t_index[i],:,loci,locj]
    qt = ql + qv
    if average_flag:
        qt = np.average(qt,axis=0)
    return qt

def get_eddy_diffusivity_flux(f,z,t_index,loci,locj,scalar,average_flag):
    '''Get eddy diffusivity flux (-K dphi/dz)
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    z: array of floats
        model height
    t_index: int, or array of ints
        index/indices for time
    loci, locj: int
        index to lat and lon
    scalar: np.array
        Scalar can be liquid water potential temperature, mixing ratio, etc
    average_flag: boolean
        whether to avearage time indicies when len(t_index)>1
    Returns
    -------
    ED: np.array
        Eddy diffusivity part of the turbulent flux [W/m^2]
    '''
    #First destagger eddy diffusivity because it's in staggered grid whereas the scalars are in regular grid
    if type(t_index) is int: #only 1 single time index is specified
        K_destagger = destagger(f['EXCH_H'][t_index,:,loci,locj],0)
        dscalar_dz = dArray_dz_2nd_central(scalar,z)
        average_flag = False #no average allowed because only 1 t_index specified
        ED = -K_destagger * dscalar_dz
    else: #t_index is a list or array
        dscalar_dz = np.zeros(len(z)-1)
        ED = np.zeros((len(t_index),len(z)))
        for tc in range(len(t_index)):
            K_destagger = destagger(f['EXCH_H'][t_index[tc],:,loci,locj],0)
            dscalar_dz = dArray_dz_2nd_central(scalar[tc,:],z)
            ED[tc,:] = -K_destagger * dscalar_dz
    if average_flag:
        ED = np.average(ED,axis=0)
    return ED, K_destagger

def dArray_dz_2nd_central(array,z):
    '''This method calculates dArray/dz using 2nd order central difference, forward/backward at end point
    Parameters
    ----------
    array: np.array
        Arrays to be differentiate
    z: np.array
        z coordinate, same size as array
    Returns
    -------
    dArray_dz: np.array
        should have the same size as array and z
    '''
    dArray_dz = np.zeros(len(z))
    for i in range(len(z)):
        if i==0: #forward difference at the start
            dArray_dz[i] = 1./np.diff(z)[i] * (array[i+1]-array[i])
        elif i==(len(z)-1): #backward difference at the end
            dArray_dz[i] = 1./np.diff(z)[i-1] * (array[i]-array[i-1])
        else: #2nd central difference for everything in the middle
            dArray_dz[i] = 1./(2.*np.diff(z)[i]) * (array[i+1]-array[i-1])   
    return dArray_dz
def get_qt(f,loci,locj):
    '''Get total water mixing ratio
    (40, 96 ocean point, 30km south of UCSD
    UCSD 53, 99)
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    loci, locj: int
        index to lat and lon
    Returns
    -------
    qt: np.array
        Total water mixing ratio
    '''
    qt = []
    ts = UCSD_WRF_timestamps2array(f)
    for i in range(1,len(ts)):
        z = getvar(f,'z',timeidx=i)[:,loci,locj].values
        t = getvar(f,'temp',timeidx=i)[:,loci,locj].values
        pblh = get_IBH_zhong_single_index(t,z)
        BL_idx = np.where(z<=pblh)[0]
        sum_of_q = f['QCLOUD'][i,:,loci,locj]+f['QRAIN'][i,:,loci,locj]+f['QVAPOR'][i,:,loci,locj]
        if len(BL_idx>0): 
            qt.append(np.average(sum_of_q[0:BL_idx[-1]+1]))
        else:
            qt.append(np.average(sum_of_q[0]))
    return np.array(qt)
def get_thetaL(f,loci,locj):
    '''Get liquid water potential temperature
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    loci, locj: int
        index to lat and lon
    Returns
    -------
    theta_l: np.array
        Liquid water potential temperature
    '''
    Lv = 2.502E6
    Rd = 287.05
    Cp = 1006.
    theta_l_ts = []
    ts = UCSD_WRF_timestamps2array(f)
    for i in range(1,len(ts)):
        p = getvar(f,'pressure',timeidx=i)[:,loci,locj].values
        z = getvar(f,'z',timeidx=i)[:,loci,locj].values
        theta = getvar(f,'theta',timeidx=i)[:,loci,locj].values
        pi = (p/p[0])**(Rd/Cp)
        ql = f['QCLOUD'][i,:,loci,locj]
        theta_l = theta - (1./pi)*(Lv/Cp)*ql
        t = getvar(f,'temp',timeidx=i)[:,loci,locj].values
        pblh = get_IBH_zhong_single_index(t,z)#f['PBLH'][i,loci,locj]
        BL_idx = np.where(z<=pblh)[0]
        print i, pblh
        if len(BL_idx>0): 
            theta_l_ts.append(np.average(theta_l[0:BL_idx[-1]+1]))
        else:
            theta_l_ts.append(np.average(theta_l[0]))
    return np.array(theta_l_ts)

def get_LWP(f,loci,locj):
    '''Get liquid water path
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    loci, locj: int
        index to lat and lon
    Returns
    -------
    LWP: np.array
        Liquid water path
    '''
    LWP = []
    qcloud = f['QCLOUD']
    ts = UCSD_WRF_timestamps2array(f)
    for i in range(1,len(ts)):
        z = getvar(f,'z',timeidx=i)[:,loci,locj]
        LWP.append(np.trapz(qcloud[i,:,loci,locj],z)*1.2*1000.)
    return np.array(LWP)

def get_cloud_top_base(f,loci,locj):
    '''Get cloud top and cloud base height
    Parameters
    ----------
    f: netCDF4.Dataset
        opened netCDF file
    loci, locj: int
        index to lat and lon
    Returns
    -------
    cloud_top, cloud_base: np.array
        Cloud top and cloud base height
    '''
    cloud_top = []
    cloud_base = []
    qcloud = f['QCLOUD']
    ts = UCSD_WRF_timestamps2array(f)
    for i in range(1,len(ts)):
        z = getvar(f,'z',timeidx=i)[:,loci,locj]
        nonzero_idx = np.nonzero(qcloud[i,:,loci,locj])[0]
        if len(nonzero_idx)>0:
            temp = z[nonzero_idx].values
            cloud_top.append(temp[temp<3000.][-1])
            cloud_base.append(z[nonzero_idx[0]])
        else:
            cloud_top.append(np.nan)
            cloud_base.append(np.nan)
    return np.array(cloud_top), np.array(cloud_base)   

def get_zone_mask(zone):
    '''
    Parameters
    ----------
    zone: str
        Zone name according to SDGE climate zones.
    Returns
    -------
    mask: np.array
    '''
    g = Dataset('/home/elynn/Documents/Ground/Climate_AvMaxTemp/WRF_d02_mask.nc')
    nc = g[zone][:,:]
    mask = ma.getmaskarray(nc)
    return mask

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

'''
These functions use octave to call matlab scripts, uncomment if needed
TMP_Inversion_Strength_Cal_V1 script in Sc-utils/Soundings
'''
# def get_IBH_zhong(t,h):
#     IBH = np.zeros((150,150))
#     for i in range(150):
#         for j in range(150):
#             hght_base = octave.TMP_Inversion_Strength_Cal_V1(t[:,i,j],h[:,i,j]/1000.,h[0,i,j]/1000.)
#             IBH[i,j] = hght_base
#     return IBH
# def get_IBH_zhong_single_index(t,h):
#     hght_base = octave.TMP_Inversion_Strength_Cal_V1(t[:],h[:]/1000.,h[0]/1000.)
#     return hght_base

def get_IBH(tc,z):
    '''Get inversion base height
    '''
    IBH = np.zeros((150,150))
    for i in range(150):
        for j in range(150):
            index = np.where(np.diff(tc[:,i,j])>0)[0][0]
            IBH[i,j] = z[index,i,j]
    return IBH
    
def find_nearest(lat,lon,value,value2):
    """Find the nearest indicies for a given location.
    
    Parameters
    --------
    lat: array_like
        latitude array
    lon: array_like
        longitude array
    value: float
        latitude for given location
    value2: float
        longitude for given location
        
    Returns
    --------
    array_like
        Indices (2D) to the nearest location.
    """
    temp_lat = lat.flatten()
    temp_lon = lon.flatten()
    idx = (np.abs(temp_lat-value)+np.abs(temp_lon-value2)).argmin()
    result = np.where(lat[:,:] == temp_lat[idx])
    return result#array[idx]
    
def UCSD_WRF_TIMESTAMP(array):
    """Process UCSD WRF timestamp to python datetime.
    
    Parameters
    --------
    array: array_like
        UCSD WRF output file variable: Times
    
    Returns
    --------
    datetime object
    """
    year = ''.join([array[0],array[1],array[2],array[3]])
    month = ''.join([array[5],array[6]])
    day = ''.join([array[8],array[9]])
    hour = ''.join([array[11],array[12]])
    min = ''.join([array[14],array[15]])
    result = datetime.datetime(int(year),int(month),int(day),int(hour),int(min))#,tzinfo=utc)#-datetime.timedelta(hours=8)
    return result

def UCSD_WRF_timestamps2array(f):
    """WRF timestamps to array using UCSD_WRF_TIMESTAMP
    Parameters
    --------
    f: opened netCDF file
        .nc file opened from Dataset
    Returns
    --------
    ts_array: array of datatime object
        Time stamps in an array.
    """    
    ts_array = []
    ts = f.variables['Times'][:,:]
    for i in range(ts.shape[0]):
        ts_array.append(UCSD_WRF_TIMESTAMP(ts[i,:]))
    ts_array = np.array(ts_array)
    return ts_array
def combine_UCSD_WRF_forecast(site,dir):
    """Combine post-processed WRF forecasts.
    
    Parameters
    --------
    site: string
        name of the site (file name)
    
    dir: string
        directory of all the post-processed csv files
        
    Returns
    --------
    None. Saved as new csv.
    
    """
    i = 0
    for file in os.listdir(dir):
        if fnmatch.fnmatch(file,site+'*.csv'):
            if i==0:
                output = pd.read_csv(dir+file,index_col=0,parse_dates=True)
            else:
                f = pd.read_csv(dir+file,index_col=0,parse_dates=True)
                output = pd.concat([output,f],axis=1)
        i = i+1
    output = output.resample('H',how='first')
    output.to_csv(dir+site+"_combined_forecast_escalator.csv")
    
    
def wrf_tk(P,theta):
    """Calculates temperature in [K] from ARW WRF model output. 
    Parameters
    --------
    P: 2D matrix
        Pressure 2D (time,height)
    theta: 2D matrix
        Potential temperature (K) 2D (time,height)
    Returns
    --------
    T: 2D matrix
        Temperature: 2D matrix (time,height)      
    """
    T = np.zeros(P.shape)
    R = 287.04
    Cp = 1004.64
    for i in range(P.shape[0]):
        T[i,:] = theta[i,:]/(P[i,0]/P[i,:])**(R/Cp)
    return T
    
def calculate_LWP_and_Inversion(f,h,location):
    """Calculates liquid water path (g/m^2) from ARW WRF model output as well as inversion strength (K) and iversion base height (km).
    Parameters
    ----------
    f: opened netCDF file
        .nc file opened from Dataset
    h: array
        height coordinate at this location
    location: tuple
        returned from the method "find_nearest"
    Returns
    -------
    LWP: array
        Liquid water path (g/m^2). The length is the same as the number of time stamps in the netCDF file
    inv_strength: array
        Inversion strength (K). 
    IBH: array
        Inversion base height (km).
    """
    qcloud = f.variables['QCLOUD'][:,:,location[0][0],location[1][0]] #cloud water mixing ratio
#    PHB = f.variables['PHB'][:,:,location[0][0],location[1][0]] #perturbation geopotential
#    PH = f.variables['PH'][:,:,location[0][0],location[1][0]] #base-state geopotential
#    h = (PH + PHB)/9.8 #from geopotential to height
#    h = (h[:,0:-1]+h[:,1:])/2. #interpolate to have the same grid points as qcloud
#    T = f.variables['T'][:,:,location[0][0],location[1][0]] #perturbation potential temperature
    P = f.variables['P'][:,:,location[0][0],location[1][0]] #perturbation pressue
    PB = f.variables['PB'][:,:,location[0][0],location[1][0]] #base-state pressure
#    T = T+300. #to go from perturbation to actual potential temperature
    P = P + PB #actual pressure
#    t = wrf_tk(P,T) #get actual temperature from pressure and potential temperature
    t = f.variables['T_PHY'][:,:,location[0][0],location[1][0]]
    qcloud = qcloud*1E3*P/(287*t) #cloud water mixing ration times density at different level
    LWP = np.zeros(qcloud.shape[0])
    inv_strength = np.zeros(qcloud.shape[0])
    IBH = np.zeros(qcloud.shape[0])
    for i in range(qcloud.shape[0]):
        LWP[i] = np.trapz(qcloud[i,:],x = h[:]) #liquid water path is the integral of vertical qcloud
#        DT_max,numofInv,hght_top,hght_base = octave.TMP_Inversion_Strength_Cal_V1(t[i,:],h[i,:]/1000.,h[i,0])
        DT_max,numofInv,hght_top,hght_base = octave.TMP_Inversion_Strength_Cal_V1(t[i,:],h[:]/1000.,h[0])
        inv_strength[i] = DT_max
        IBH[i] = hght_base
    return LWP, inv_strength, IBH
    

def calculate_zone_LWP_and_Inversion(f,mask):
    """Calculates liquid water path (g/m^2) from ARW WRF model output as well as inversion strength (K) and iversion base height (km).
    Parameters
    ----------
    f: opened netCDF file
        .nc file opened from Dataset
    h: array
        height coordinate at this location
    mask: np.array
        mask of true/false for the given climate zone
    Returns
    -------
    LWP: array
        Liquid water path (g/m^2). The length is the same as the number of time stamps in the netCDF file
    inv_strength: array
        Inversion strength (K). 
    IBH: array
        Inversion base height (km).
    """
    qcloud = f.variables['QCLOUD']
    qcloud = select_zone_by_mask(qcloud,mask)
    P = f.variables['P'] #perturbation pressue
    P = select_zone_by_mask(P,mask)
    PB = f.variables['PB'] #base-state pressure
    PB = select_zone_by_mask(PB,mask)
    PHB = f.variables['PHB'] #perturbation geopotential
    PHB = select_zone_by_mask(PHB,mask)
    PH = f.variables['PH'] #base-state geopotential
    PH = select_zone_by_mask(PH,mask)
    h = (PH + PHB)/9.8 #from geopotential to height
    h = (h[0,0:-1]+h[0,1:])/2. #interpolate to have the same grid points as qcloud
    P = P + PB #actual pressure
    t = f.variables['T_PHY']
    t = select_zone_by_mask(t,mask)
    qcloud = qcloud*1E3*P/(287*t) #cloud water mixing ration times density at different level
    LWP = np.zeros(qcloud.shape[0])
    inv_strength = np.zeros(qcloud.shape[0])
    IBH = np.zeros(qcloud.shape[0])
    for i in range(qcloud.shape[0]):
        LWP[i] = np.trapz(qcloud[i,:],x = h[:]) #liquid water path is the integral of vertical qcloud
        DT_max,numofInv,hght_top,hght_base = octave.TMP_Inversion_Strength_Cal_V1(t[i,:],h[:]/1000.,h[0])
        inv_strength[i] = DT_max
        IBH[i] = hght_base
    return LWP, inv_strength, IBH

def select_zone_by_mask(var,mask):
    '''This method takes in a variable and a 2d zone mask, then average all the points inside the zone.
    Parameters
    ----------
    var: opened netCDF variable
        variable opened from .nc file
    mask: np.array
        2d array of true/false points (true for points inside given climate zone)
    Returns
    -------
    var: np.array
        zone averaged variable
    '''
    ndim = len(var.shape) #number of dimensions
    if ndim == 3:
        mask = np.repeat(mask[np.newaxis,:,:],var.shape[0],axis=0)
        var = ma.masked_where(mask,var)
        var = ma.average(var,axis=2)
        var = ma.average(var,axis=1)        
    elif ndim == 4:
        mask = np.repeat(mask[np.newaxis,:,:],var.shape[1],axis=0)
        mask = np.repeat(mask[np.newaxis,:,:],var.shape[0],axis=0)
        var = ma.masked_where(mask,var)
        var = ma.average(var,axis=3)
        var = ma.average(var,axis=2)
    return var
    
def get_WRF_3d_variable(f,location,var_string):
    """Get 3D WRF variable.
    Parameters
    ----------
    f: opened netCDF file
        .nc file opened from Dataset
    location: tuple
        returned from the method "find_nearest"
    var_string: str
        Variable name (e.g. T2 for temp at 2m, SWDOWN for GHI)
    Returns
    -------
    var: array
        Variable. The length is the same as the number of time stamps in the netCDF file.
    """
    var = f.variables[var_string][:,location[0][0],location[1][0]]
    return var

def get_WRF_zone_3d_variable(f,mask,var_string):
    """Get 3D WRF variable.
    Parameters
    ----------
    f: opened netCDF file
        .nc file opened from Dataset
    mask: np.array
        mask of true/false for the given climate zone
    var_string: str
        Variable name (e.g. T2 for temp at 2m, SWDOWN for GHI)
    Returns
    -------
    var: array
        Variable. The length is the same as the number of time stamps in the netCDF file.
    """
    var = f.variables[var_string]
    var = select_zone_by_mask(var,mask)
    return var
    
def WRF_to_30min_ending(wrf):
    """Shift WRF data to 30min ending to match observation.
    Parameters
    ----------
    wrf: pandas dataframe
        Input WRF data (15min instantaneous)
    Returns
    -------
    wrf: pandas dataframe
        Output 30min ending WRF
    """
    wrf2 = wrf.resample('1min').mean().interpolate(limit=15)
    wrf2 = wrf.shift(-15).copy()
    wrf2 = wrf2.resample('30min',label='right',closed='right').mean()
    return wrf

def UTC_to_PPT(wrf):
    """Shift WRF data from UTC to PPT (with daylight savings)
    Parameters
    ----------
    wrf: pandas dataframe
        Input WRF data (time resolution: freq)
    freq: int
        Time resolution
    Returns
    -------
    new_wrf: pandas dataframe
        Output WRF data in PPT
    """
    ts = wrf.index
    ts_utc = ts.tz_localize(pytz.UTC)
    ts_pacific = []
    pacific = pytz.timezone('US/Pacific')
    for i in range(len(ts)):
        ts_pacific.append(ts_utc[i].astimezone(pacific).tz_localize(None))
    ts_pacific = np.array(ts_pacific)
    new_wrf = pd.DataFrame(wrf.as_matrix(),index=ts_pacific,columns=wrf.columns)
    new_wrf = new_wrf.resample('H').first()
    return new_wrf
