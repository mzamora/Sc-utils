#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Stitch the raw files from ucla-les from different processors
This script has only been tested on nx = 2 and ny = 4, use with caution
@author: elynn
"""

import numpy as np
from netCDF4 import Dataset
'''Configure here'''
file_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/'
output_prefix = 'rf01'
output_dir = '/home/elynn/Documents/uclales/dycomsrf01_br02/RF01_hour34_1min/'
nx = 2 #number of processors in x-direction
ny = 4 #number of processors in y-direction
x_pts = 48 #number of points in x-direciton for each processor
y_pts = 24 #number of points in y-direciton for each processor
z_pts = 131 #number of points in z-direction
t_pts = 60 #number of time index 
variables = ['l','w','t','q','u','v'] #variables to be stitched
'''End configuration'''

domain_var = {} #dictionary to store all variables
for i in range(len(variables)): #initialize the stitched variables
    domain_var[variables[i]] = np.zeros((t_pts,y_pts*ny,x_pts*nx,z_pts))
xm_all = []
ym_all = []
for x in range(nx):
    for y in range(ny):
        filename = file_dir+output_prefix+'.'+str(x).zfill(4)+str(y).zfill(4)+'.nc' #filename from each processor
        current = Dataset(filename)
        for i in range(len(variables)):
            domain_var[variables[i]][:,y*y_pts:y*y_pts+y_pts,x*x_pts:x*x_pts+x_pts,:] = current[variables[i]][:,:,:,:]
        zt = current['zm'][:]
        xm_all.append(current['xm'][:])
        ym_all.append(current['ym'][:])
        current.close()
        
x = sorted(np.unique(np.array(xm_all))) #get x grids from processors
y = sorted(np.unique(np.array(ym_all))) #get y grids from processors
time = np.arange(1,t_pts+1)

root_grp = Dataset(output_dir+output_prefix+'_stitched1.nc', 'w')
root_grp.description = 'Stitched from 8 different processors'

# dimensions
root_grp.createDimension('time', t_pts)
root_grp.createDimension('xm', nx*x_pts)
root_grp.createDimension('ym', ny*y_pts)
root_grp.createDimension('zm', z_pts)

# variables
times = root_grp.createVariable('time', 'f4', ('time'))
xcoord = root_grp.createVariable('xm', 'f4', ('xm'))
ycoord = root_grp.createVariable('ym', 'f4', ('xm'))
zcoord = root_grp.createVariable('zm', 'f4', ('zm'))
temps = {} #initialize netCDF4 variables to be written
for i in range(len(variables)):
    temps['temp'+str(i+1)] = root_grp.createVariable(variables[i], 'f4', ('time','ym', 'xm', 'zm'))
    temps['temp'+str(i+1)][:,:,:,:] = domain_var[variables[i]][:,:,:,:]
times[:] = time[:]
xcoord[:] = x[:]
ycoord[:] = y[:]
zcoord[:] = zt[:]
root_grp.close()