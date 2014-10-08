# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/eastwoodknight/.spyder2/.temp.py
"""
import pandas as pd
from fit_exp import fit_exp

data = []
frame1 = pd.read_table('data2.txt',skiprows=1,header=None,nrows=10)
frame2 = pd.read_table('data2.txt',skiprows=12,header=None,nrows=10)
frame3 = pd.read_table('data2.txt',skiprows=23,header=None,nrows=10)
data.append(frame1)
data.append(frame2)
data.append(frame3)

def func(p1,x1):
    return p1[0]*np.ones(len(x1)) + p1[1]*np.exp(-x1/p1[2])
  
p0=(2100,1350,-1110)
  
f = open('Clbr_coefficients','w')
f.write('%-7s\t%-7s\t%-7s\t%-7s\t\n' %('№','A','B','C'))
number = 0
for frame in data:
    for i in frame1:
        number += 1
        if i > 0:
            x = frame[i]
            en = frame[0]
            coef = fit_exp(p0,x,en,func,show=False)
            s = '%-5d\t%-5.1f\t%-5.1f\t%-5.1f\t\n'%(number,coef[0],coef[1],coef[2])
            f.write(s)
f.close()
            

