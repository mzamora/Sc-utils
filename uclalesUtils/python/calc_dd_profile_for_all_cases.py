#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 18:07:11 2018

@author: elynn
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)

'''
For downdraft theta difference
'''
alpha = 0.98
mu_w = 0.59
sigma_wi = 0.36
vertical_fluxes = pd.read_csv('/home/elynn/Documents/uclales/dycomsrf01_br02_no_mean_uv/hr_3_4_1min/Vertical_fluxes_all_time.csv',index_col=0)
vertical_fluxes = vertical_fluxes[vertical_fluxes.columns[0::4]].mean(axis=1)
PBL_flux = vertical_fluxes.loc[0.96:1].mean()
theta_diff = alpha*PBL_flux/sigma_wi
print 'DYCOMS RF01 no wind:', theta_diff, 'K'

vertical_fluxes = pd.read_csv('/home/elynn/Documents/uclales/cgils_s12_ctl/hr_3_4_1min/plume_output/zi2_bar_5percent_plume/Vertical_fluxes_all_time.csv',index_col=0)
vertical_fluxes = vertical_fluxes[vertical_fluxes.columns[0::4]].mean(axis=1)
PBL_flux = vertical_fluxes.loc[0.96:1].mean()
theta_diff = alpha*PBL_flux/sigma_wi
print 'CGILS s12:', theta_diff, 'K'

'''
For downdraft qt difference
'''
alpha = -1.817
mu_w = 0.59
sigma_wi = 0.36
vertical_fluxes = pd.read_csv('/home/elynn/Documents/uclales/dycomsrf01_br02_no_mean_uv/hr_3_4_1min/Vertical_fluxes_all_time.csv',index_col=0)
vertical_fluxes = vertical_fluxes[vertical_fluxes.columns[2::4]].mean(axis=1)
PBL_flux = vertical_fluxes.loc[0.96:1].mean()
qt_diff = alpha*PBL_flux/sigma_wi
print 'DYCOMS RF01 no wind:', qt_diff*1000., 'g/kg'

vertical_fluxes = pd.read_csv('/home/elynn/Documents/uclales/cgils_s12_ctl/hr_3_4_1min/plume_output/zi2_bar_5percent_plume/Vertical_fluxes_all_time.csv',index_col=0)
vertical_fluxes = vertical_fluxes[vertical_fluxes.columns[2::4]].mean(axis=1)
PBL_flux = vertical_fluxes.loc[0.96:1].mean()
qt_diff = alpha*PBL_flux/sigma_wi
print 'CGILS s12:', qt_diff*1000., 'g/kg'