#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 16:08:42 2018

@author: elynn
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import Sc_methods as Sc_m
from netCDF4 import Dataset
import seaborn as sns #just to make plots look nicer, optional
sns.set(style='ticks',font_scale=1.2)
from oct2py import octave
octave.addpath('/home/elynn/Documents/github/Sc-utils/Soundings')

z = np.linspace(0,1101,10)
sigma_w1 = Sc_m.holtslag_Moeng_91(1.5,0.07,900.,z)
sigma_w2 = Sc_m.holtslag_Moeng_91(1.4,0.02,900.,z)
sigma_w3 = Sc_m.holtslag_Moeng_91(1.7,0.03,500.,z)
plt.plot(sigma_w1,z/900.,label=r'$\sigma_{w1}$')
plt.plot(sigma_w2,z/900.,label=r'$\sigma_{w2}$')
plt.plot(sigma_w3,z/500.,label=r'$\sigma_{w3}$')
plt.legend(loc=0)