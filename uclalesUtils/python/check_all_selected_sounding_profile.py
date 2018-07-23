#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 11:04:14 2018

@author: elynn
"""
import sys
sys.path.append('/mnt/lab_45d1/database/Sc_group')
import numpy as np
import prep_uclales_methods as prep_les
import matplotlib.pyplot as plt
import pandas as pd
dates = pd.read_csv('/mnt/lab_45d1/database/Sc_group/NKX_sounding/Final_testing_dates_selection.csv',index_col=0,parse_dates=True)
home_dir = '/mnt/lab_45d1/database/Sc_group/'
for i in range(dates.shape[0]):
    year = str(dates.index[i].year)
    month = str(dates.index[i].month)
    day = str(dates.index[i].day)
    date = pd.date_range(year+'-'+month+'-'+day+' 12:00', year+'-'+month+'-'+day+' 12:00')[0]
    caseName = 'NKX_'+date.strftime('%Y%m%d')
    IC, uwind, vwind, z = prep_les.make_sound_in_file(home_dir,date)
    f = np.genfromtxt('sound_in',skip_header=1)
    plt.plot(f[:,1],f[:,0],'-o')
    plt.title(year+'-'+month+'-'+day)
    plt.savefig('/mnt/lab_45d1/database/Sc_group/LES_analysis/initial_sounding_check/'+date.strftime('%Y%m%d')+'_idealized_temperature_profile.png',dpi=200,bbox_inches='tight')
    plt.close()