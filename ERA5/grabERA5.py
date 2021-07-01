#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grabbing data from ERA 5 for Arica, Chile
M. Zamora Zapata - U. Chile
mzamora.github.io
"""

# Grabbing data from ERA 5
# ref: https://towardsdatascience.com/read-era5-directly-into-memory-with-python-511a2740bba0

import cdsapi
import xarray as xr
from urllib.request import urlopen
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import calendar

c = cdsapi.Client()

# dataset to read: vertical profiles
dataset = 'reanalysis-era5-pressure-levels'

# flag to download data
download_flag = False

# location: arica's airport
lat0=-18.3514
lon0=-70.3358
ele0=59
preslevels=[1, 2, 3, 5, 7, 10, 20, 30, 50, 70]
preslevels.extend(list(range(100, 1000, 25)))

# set dates
year=2017
ntime=len(pd.date_range(pd.datetime(year,1,1),pd.datetime(year,12,31,23,0,0),freq='H')) #length of expected time object for the year
nz=36 #len(preslevels) #expected height variable length
it=0

# initialize empty matrices to fill
qtout=np.empty((ntime,nz))
tempout=np.empty((ntime,nz))
uout=np.empty((ntime,nz))
vout=np.empty((ntime,nz))
divout=np.empty((ntime,nz))
qlout=np.empty((ntime,nz))
rhout=np.empty((ntime,nz))
ccout=np.empty((ntime,nz))
timeout=np.empty(ntime,dtype='datetime64[ns]')

for month in range(1,13): # go through all year, day by day
    date0=pd.datetime(year,month,1)
    date1=pd.datetime(year,month,calendar.monthrange(year,month)[1])
    wanteddates=pd.date_range(date0,date1)
    for datei in wanteddates:
        # api parameters 
        params = {
            'format': 'netcdf',
            'product_type': 'reanalysis',
            'variable': ['specific_humidity','temperature','u_component_of_wind',
                         'v_component_of_wind','vertical_velocity','divergence',
                         'relative_humidity','specific_cloud_liquid_water_content',
                         'fraction_of_cloud_cover'],
            'pressure_level':preslevels,
            'year': [datei.strftime('%Y')], #['2020'],
            'month': [datei.strftime('%m')], #['01'], [str(i).zfill(2) for i in range(0,13)], 
            'day':  [datei.strftime('%d')], #['01'],[str(i).zfill(2) for i in range(0,32)],
            'time': [str(i).zfill(2)+':00' for i in range(0,24)],
           'grid': [0.25, 0.25],
            'area': [lat0-.1, lon0-.1, lat0+.1, lon0+.1], #.1 gives single point data
        }
    
        # retrieves the path to the file
        fl = c.retrieve(dataset, params)
        
        # download the file 
        if download_flag:
            fl.download("./output.nc")
            
        # load into memory
        with urlopen(fl.location) as f:
            ds = xr.open_dataset(f.read())
        
        t=ds['time'].data
        nt=len(t)
        
        timeout[it:it+nt]=t
        qtout[it:it+nt,:]=ds['q'].data.squeeze()
        tempout[it:it+nt,:]=ds['t'].data.squeeze()
        uout[it:it+nt,:]=ds['u'].data.squeeze()
        vout[it:it+nt,:]=ds['v'].data.squeeze()
        divout[it:it+nt,:]=ds['d'].data.squeeze()
        rhout[it:it+nt,:]=ds['r'].data.squeeze()
        qlout[it:it+nt,:]=ds['clwc'].data.squeeze()
        ccout[it:it+nt,:]=ds['cc'].data.squeeze()
        
        it=it+nt
        print('===========')
        print(datei)
        print('===========')
    
    # store results every month for safety
    #timeout.tofile('ERA5/time_'+str(year)+'.csv',sep=',')
    np.array(timeout).tofile('ERA5/time_'+str(year)+'.csv',sep=',')
    qtout.tofile('ERA5/qt_'+str(year)+'.csv',sep=',')
    tempout.tofile('ERA5/temp_'+str(year)+'.csv',sep=',')
    uout.tofile('ERA5/u_'+str(year)+'.csv',sep=',')
    vout.tofile('ERA5/v_'+str(year)+'.csv',sep=',')
    divout.tofile('ERA5/div_'+str(year)+'.csv',sep=',')
    rhout.tofile('ERA5/rh_'+str(year)+'.csv',sep=',')
    qlout.tofile('ERA5/ql_'+str(year)+'.csv',sep=',')
    ccout.tofile('ERA5/cc_'+str(year)+'.csv',sep=',')

pres=ds['level'].data.squeeze()
pres.tofile('ERA5/preslevels.csv',sep=',')



