# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 14:06:39 2014

@author: eastwood
"""
import pandas as pd
from fit_exp import fit_exp
    
def read_coefs(filename):
    table = pd.read_table(filename,sep='  ',header=None,nrows = 48)
    f_parse = lambda x: float(x.split()[1] )
    coef_a = table[1].map(f_parse)
    coef_b = table[2].map(f_parse)
    coef_a.index= np.arange(1,49)
    coef_b.index= np.arange(1,49)
    frame_front = pd.DataFrame()
    frame_front['A'] = coef_a
    frame_front['B'] = coef_b
    
    table = pd.read_table(filename,sep='  ',header=None,skiprows=49,nrows = 128)
    coef_a = table[1].map(f_parse)
    coef_b = table[2].map(f_parse)
    coef_a.index= np.arange(1,129)
    coef_b.index= np.arange(1,129)
    frame_back = pd.DataFrame()
    frame_back['A'] = coef_a
    frame_back['B'] = coef_b
    
    return frame_front, frame_back
    
def func(p1,x1):
    return p1[0]*np.ones(len(x1)) + p1[1]*np.exp(-x1/p1[2])
    
def process(channel,energy,strip,coefs):
    #calculate energies of last 5 numbers in Y column
    energy1 = np.array(energy)
    channel1 = np.array(channel)
    energy1[8] = channel1[8]*coefs['A'][strip] + coefs['B'][strip]
    energy1[9] = energy1[8]*2
    energy1[10] = energy1[8]*10
    energy1[11] = energy1[8]*20
    energy1[12] = energy1[8]*28
    #fit by exp
    coef = fit_exp(p0,channel1,energy1,func,show=False)
    return coef,channel1,energy1
				
def process1(channel,energy,strip):
    energy1 = np.array(energy)[:7]
    channel1 = np.array(channel)[:7]
    p0=(6000,450,-1000)
    coef1 = fit_exp(p0,channel1,energy1,func,show=False)			
    energy2 = np.array(energy)[7:]
    channel2 = np.array(channel)[7:]
    p0 = (8000,-8000,-100)			
    coef2 = fit_exp(p0,channel2,energy2,func,show=False)
    return coef1,channel1,energy1,	coef2,channel2,energy2
				
def process2(channel,energy,strip):
    #print len(channel),len(energy)
    energy1 = np.array(energy)
    channel1 = np.array(channel)
    coef1 = fit_exp(p0,channel1,energy1,func,show=False)			
    
    return coef1,channel1,energy1

#def read_table(filename): 
#    return pd.read_table(filename,sep='\s+',skiprows=1)
    """
    data = []
    data.append( pd.read_table(filename,skiprows=1,nrows=13) )
    data.append( pd.read_table(filename,skiprows=14,nrows=13) )
    data.append( pd.read_table(filename,skiprows=27,nrows=13) )
    return data
    """
    
def write_strip_report(Fobj,strip,coef,channel,energy):
    coef = coef[0] 
    Fobj.write('Strip%-2s %d %2s= %4.2f ; %2s= %4.2f ; %2s= %4.2f \n' %('№',strip,'A',coef[0],'B',coef[1],'C',coef[2]))  
    Fobj.write('  %-5s %-6s   %-11s' %('chan','energy','calc_energy \n'))
    energy1 = func(coef,channel)    
    for i in xrange(len(channel)):
         Fobj.write('  %-4d  %-6d   %-6.2f \n' %(channel[i],energy[i],energy1[i]) )
    Fobj.write('\n')
    
if __name__ == '__main__':
    
    p0=(6000,450,-1000)
    #p0 = (8000,-8000,-100)
    
    list_front = ['Sheet1.t','Sheet2.t','Sheet3.t']
    #list_back = ['Sheet4.t','Sheet5.t','Sheet6.t','Sheet7.t','Sheet8.t','Sheet9.t','Sheet10.t','Sheet11.t'] 
    list_back = ['Sheet1b.t','Sheet2b.t','Sheet3b.t','Sheet4b.t','Sheet5b.t','Sheet6b.t','Sheet7b.t','Sheet8b.t']
    f = open('Calibration_report_24Oct2014.txt','w')

    A,B,C = [],[],[]
    f.write('Front detectors coefficients\n')
    ind = np.arange(1,17)
    for filename in list_front:
        frame = pd.read_table(filename,sep='\s+',skiprows=1)
        for i,col in enumerate(frame.columns[1:]):
			indx = frame[col].notnull()
			if indx.sum() == 0:
				continue
			else:
				coef, channel, energy = process2(frame[col][indx],frame['E'][indx],ind[i])
				A.append(coef[0][0])
				B.append(coef[0][1])
				C.append(coef[0][2])
				#print '\n coef:',coef, channel, energy, coef2,channel2,energy2
				write_strip_report(f,ind[i],coef,channel,energy)
				#write_strip_report(f,ind[i],coef2,channel2,energy2)
				#f.write('\n\n')
        ind += 16  
    

    f.write('\n\n')				
    for Lst in [A,B,C]:
        f.write('\n\n')
        k = 0
        for j in Lst:
            f.write('%4.2f,' %(j) )
            k += 1
            if k == 6:
			   f.write('\n')
			   k = 0
    f.close()    
	
							
#"""			
				
"""	
    f.write('Back detectors coefficients\n')
    ind = np.arange(1,17)
    for filename in list_back:
        frame = pd.read_table(filename,sep='\s+',skiprows=1)
        for i,col in enumerate(frame.columns[1:]):
			coef, channel, energy= process2(frame[col],frame['E'],ind[i])
			#print '\n coef:',coef, channel, energy, coef2,channel2,energy2
			write_strip_report(f,ind[i],coef,channel,energy)
        ind += 16    
    f.close()
"""
""" September version
	coefs_front,coefs_back = read_coefs('AB.txt')
	
    f.write('Front detectors coefficients\n')
    ind = np.arange(1,17)
    for filename in list_front:
        frame = pd.read_table(filename,sep='\s+',skiprows=1)
        for i,col in enumerate(frame.columns[1:]):
            coef, channel, energy = process(frame[col],frame['E'],ind[i],coefs_front)
            write_strip_report(f,ind[i],coef,channel,energy)
        ind += 16
        
    f.write('***************************\n')
    f.write('Back detectors coefficients\n')
    ind = np.arange(1,17)
    for filename in list_back:
        frame = pd.read_table(filename,sep='\s+',skiprows=1)
        for i,col in enumerate(frame.columns[1:]):
            coef, channel, energy = process(frame[col],frame['E'],ind[i],coefs_back)
            write_strip_report(f,ind[i],coef,channel,energy)
        ind += 16
    
    f.close()
"""
