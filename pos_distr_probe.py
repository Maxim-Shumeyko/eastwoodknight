# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:19:15 2014

@author: eastwoodknight
"""

import read_files
import numpy as np
from matplotlib import pyplot as plt

def probe(flag = False):
    frame = read_files.read_files('tsn.456-tsn.461',strip_convert=flag)
    mask = (frame['tof']==0)
    hist_f = np.histogram(frame['strip'][mask],bins = np.arange(0,49))
    hist_b1 = np.histogram(frame['b_strip 1'][mask],bins=np.arange(0,129))
    hist_b2 = np.histogram(frame['b_strip 2'][mask],bins=np.arange(0,129))
    
    fig, ax = plt.subplots(2,2,sharex = True)
    ax[0,0].plot(hist_f[0],linestyle='steps')
    ax[0,1].plot(hist_b1[0],linestyle='steps')
    ax[1,1].plot(hist_b2[0],linestyle='steps')
        
    plt.show()