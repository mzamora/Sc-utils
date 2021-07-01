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

c = cdsapi.Client()

# dataset to read: surface data
dataset = 'reanalysis-era5-single-levels'

# flag to download data
download_flag = False

# location: arica's airport
lat0=-18.3514
lon0=-70.3358
ele0=59

df = pd.DataFrame([])

year=2021 #set year to download
#for date in wanteddates:
# api parameters 
params = {
'format': 'netcdf',
'product_type': 'reanalysis',
'variable': ['cloud_base_height','boundary_layer_height','total_cloud_cover',
             'low_cloud_cover','medium_cloud_cover','high_cloud_cover',
             'surface_sensible_heat_flux','surface_latent_heat_flux',
             'surface_solar_radiation_downwards','mean_sea_level_pressure'],
# we can grab the whole year in this case (less data)
'year': [str(year)], #[datei.strftime('%Y')], #['2020'], 
'month': [str(i).zfill(2) for i in range(0,5)], #[datei.strftime('%m')], #['01'],
'day': [str(i).zfill(2) for i in range(0,32)], #[datei.strftime('%d')], #['01'],
'time': [str(i).zfill(2)+':00' for i in range(0,24)],
'grid': [0.25, 0.25],
'area': [lat0-.1, lon0-.1, lat0+.1, lon0+.1],
}

# retrieves the path to the file
fl = c.retrieve(dataset, params)

# download the file
if download_flag:
    fl.download("./output.nc")

# load into memory
with urlopen(fl.location) as f:
    ds = xr.open_dataset(f.read())

t=ds['time'].data.squeeze()
zb=ds['cbh'].data.squeeze()
zi=ds['blh'].data.squeeze()
tcc=ds['tcc'].data.squeeze()
lcc=ds['lcc'].data.squeeze()
mcc=ds['mcc'].data.squeeze()
hcc=ds['hcc'].data.squeeze()
shf=ds['sshf'].data.squeeze()
lhf=ds['slhf'].data.squeeze()
sw=ds['ssrd'].data.squeeze()
slp=ds['msl'].data.squeeze()
#output dictionary
data={'zb':zb,'zi':zi, 'cc':tcc, 'lcc':lcc,'mcc':mcc, 'hcc':hcc, 'shf':shf,
      'lhf':lhf, 'sw':sw, 'slp':slp}

#weird stuff with 2 column split data in 2021 (exper variable)
t=ds['time'].data.squeeze()
zb=ds['cbh'][:,0,0,0].data.squeeze()
zi=ds['blh'][:,0,0,0].data.squeeze()
tcc=ds['tcc'][:,0,0,0].data.squeeze()
lcc=ds['lcc'][:,0,0,0].data.squeeze()
mcc=ds['mcc'][:,0,0,0].data.squeeze()
hcc=ds['hcc'][:,0,0,0].data.squeeze()
shf=ds['sshf'][:,0,0,0].data.squeeze()
lhf=ds['slhf'][:,0,0,0].data.squeeze()
sw=ds['ssrd'][:,0,0,0].data.squeeze()
slp=ds['msl'][:,0,0,0].data.squeeze()
data={'zb':zb,'zi':zi, 'cc':tcc, 'lcc':lcc,'mcc':mcc, 'hcc':hcc, 'shf':shf,
      'lhf':lhf, 'sw':sw, 'slp':slp}

df=pd.DataFrame(data,index=t)

df.to_csv('ERA5/SRF_'+str(year)+'.csv')




