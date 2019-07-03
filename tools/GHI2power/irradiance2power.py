#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 13:48:28 2018

@author: elynn
"""
import numpy as np
import datetime
import pandas as pd
import pvlib
from pvlib import pvsystem
from pvlib.location import Location
from pvlib import clearsky
from pvlib import irradiance
from pvlib import atmosphere
from pvlib.location import Location
from pvlib.modelchain import ModelChain
import seaborn as sns
sns.set(style="ticks",font_scale=1.1)

obs = pd.read_csv('Hopkins_Park_PV_2018.csv',index_col=0,parse_dates=True)
index = []
for i in range(len(obs.index)):
    index.append(pd.Timestamp(obs.index[i]+datetime.timedelta(hours=7),tz='Etc/UTC'))
obs = pd.DataFrame(obs['Real Power Mean'].values,index=index,columns=['Obs power'])
'''
load modules and inverters
'''
sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
cec_inverters = pvlib.pvsystem.retrieve_sam('cecinverter')
inverter = cec_inverters['Xantrex_Technology__GT3_8_NA_240_208_UL_05__208V__208V__CEC_2018_']
module = sandia_modules['Kyocera_Solar_KD205GX_LP__2008__E__']
system = pvsystem.PVSystem(surface_tilt=10, surface_azimuth=180,
                  module_parameters=module,
                  inverter_parameters=inverter) #only takes one single module
system.module_parameters['pdc0'] = 10
system.module_parameters['gamma_pdc'] = -0.004 #temperature efficiency loss
#1652 panels
#print system.module_parameters
sd = Location(32.883729,-117.239341)
'''input own weather data'''
weather_raw = pd.read_csv('/Users/elynn/Documents/powerForecast/WRF_intraday_weather_forecast.csv',index_col=0,parse_dates=True)
weather_raw['wind_speed'] = np.sqrt(weather_raw['U10']**2 + weather_raw['V10']**2)
weather_raw = weather_raw[['SWDOWN','SWDDNI','SWDDIF','T2','wind_speed']]
weather_raw.columns = ['ghi', 'dni', 'dhi', 'temp_air', 'wind_speed']
weather_raw['temp_air'] = weather_raw['temp_air'] - 273.15
index = []
for i in range(len(weather_raw.index)):
    index.append(pd.Timestamp(weather_raw.index[i],tz='Etc/UTC'))
weather = pd.DataFrame(weather_raw.values,index=index,columns=['ghi', 'dni', 'dhi', 'temp_air', 'wind_speed'])
mc = ModelChain(system, sd, aoi_model='physical',name='test')
mc.transposition_model = 'perez'
mc.run_model(times=weather.index, weather=weather)
power = mc.ac * 1652 /1000.
power.name = 'Power forecast'
power = power.to_frame().resample('15min').first()
'''
simulation complete
'''
obs = obs.resample('15min').first()
comp = pd.concat((obs,power),axis=1)
comp.plot()
