# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 13:07:54 2014

@author: eastwoodknight
"""

from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import numpy as np

# calibration function
"""
#    Fit by function F(x) = A + B*exp(-x/C)
"""
def fit_exp(p0,x,y,func,show=True):

    def residuals(p1,y1,x1):
        return y1 - func(p1,x1)
    """    
    def Jac(p0,y,x):
        dfdA = np.ones(len(x)) #* residuals(p0,y,x) 
        dfdB = 2*x#* residuals(p0,y,x) 
        #dfdC =-  p0[1]*dfdB/p0[2]
        deriv = np.c_[dfdA, dfdB]#,dfdC]
        return  deriv
    """
        
    solution = leastsq(residuals, p0,args = (y,x) )
    
    if show:
        fig = plt.figure()
        subp = fig.add_subplot(111)
        subp.plot(x,y,'hg',label='True data') 
        xmin = x.min()
        xmax = x.max()
        x = np.linspace(xmin,xmax,500)
        subp.plot(x,func(solution[0],x),'--r',label='Fitting function') 
        subp.legend(loc='upper left')
        plt.show()
       
    return solution
    
if __name__ == '__main__':
    def func(p1,x):
        return p1[0]*np.ones(len(x)) + p1[1]*np.exp(x/p1[2])
        
    x = np.linspace(0,50,100)
    p0 = [0, 2, 5]   
    k = 0
    y = func(p0,x) + k*np.random.rand(len(x)) - float(k)/2    
    
    coef = fit_exp(p0,x,y,func,True)   
    print coef
    
    
    