# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 17:10:17 2013

@author: eastwoodknight
"""
import pandas as pd
from pandas import Series, DataFrame
from datetime import datetime
import numpy as np
from scipy.linalg import lstsq
import matplotlib.pyplot as plt

def fit_method(x,y,show = False):
    x = np.array(x)
    #x = x.astype(np.float32)
    A = np.c_[np.ones(len(x)),x]
    solution = lstsq(A,y)
    if show:
        plt.plot(x,solution[0][1]*x + solution[0][0]*np.ones(len(x)),'r--',
                 x,y,'bD')
        plt.legend(['fit line','data points'],loc = 'upper left')
        plt.show()
    return solution[0]
        
def calibrate(filename):
    energies = [9261,8700,7923,7137,6999]
    msv_len = len(energies)
    #read data
    lines = open(filename,'r').readlines()
    clbr_data = []
    for i in range(1,len(lines)):
        line = lines[i].split()
        channels = []
        for j in range(1,1+msv_len):
            channels.append(line[j])
        clbr_data.append((channels,energies))
    #calibration
    solutions = []
    for x,y in clbr_data:
        solutions.append( fit_method(x,y,False) ) 
    #write results to file
    f = open( 'coefficients'+str(datetime.now()),'w' )
    for i in solutions:
        f.write('%.4f  %.4f \n'%(i[1],i[0])) #str(i[0])+' '+str(i[1])+'\n')
    f.close()
    return solutions
    
def calibrate_pd(filename):
    energies = [9261,8700,7923,7137]
    #read data
    df =  pd.read_table(filename,skiprows = [0,],sep = '\s+', header = None) 
    #calibration
    solutions = []
    for i in range(len(df)):
        solutions.append( fit_method(df.ix[i][1:],energies,False) ) 
    #write results to file
    f = open( 'coefficients'+str(datetime.now()).split('.')[0],'w' )
    for i in solutions:
        f.write('%.4  %.4 \n'%[i[1],i[0]])
    f.close()
    return solutions



if __name__ == '__main__':
    data = calibrate('peaks_Jan23.txt')
    