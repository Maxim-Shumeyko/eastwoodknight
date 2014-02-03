# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from pandas import Series as Ser, DataFrame as DF
import Visualisator
import time as t
"""  
######################################################################
The module is for data processing on experiment of superheavie synthesis (GFS).
It includes some functions to work with data files of binary format and
provides essential functionality to read them and to build some useful distributions.
######################################################################
FUNCTIONS:
    
1.  read_file(filename) - to read single file.
    return (pandas) DataFrame
    
2.  read_files(filenames) - to read group of files
    # filenames can be put in like this: "tsy.345-tsy.349", 'tsy.322 - tsy.330', 'tsy.322,tsy.323,tsy.324'
    return (pandas) DataFrame
    
3.  read_times(filename) - to read only time and strip information about event
    return (pandas) DataFrame
    
4.  get_front_spectrs(filename) - to get alpha front spectrs
    # data - filename or dataframe
    return list of spectrs (48 spectrs, 8192 channels each)
    - show histograms of distributions
    
5.  tof_distr(filename) - to get information about tof-distributions
    output text information
    
"""

def read_file(filename,strlength = 14,write_file=False,energy_scale=False,time_corr=False,strip_convert=True):
     data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,strlength)#[:,0:9]
     #beam marker
     beam_marker = (data[:,0]>>4)% 16
     #front detectors
     indexes = data[:,0] % 16 
     strips = np.array(data[:,2] / 4096 + 1,dtype='uint16')
     strips = np.where( indexes < 3, indexes*16+strips, strips)
     channel = data[:,1] % 8192
     
     #front strips convert func
     fs_convert = {}
     fs_convert[1]=1
     fs_convert[2]=17
     fs_convert[3]=33
     fs_convert[4]=2
     fs_convert[5]=18
     fs_convert[6]=34
     fs_convert[7]=3
     fs_convert[8]=19
     fs_convert[9]=35
     fs_convert[10]=4
     fs_convert[11]=20
     fs_convert[12]=36
     fs_convert[13]=5
     fs_convert[14]=21
     fs_convert[15]=37
     fs_convert[16]=16
     fs_convert[17]=6
     fs_convert[18]=22
     fs_convert[19]=38
     fs_convert[20]=7
     fs_convert[21]=23
     fs_convert[22]=39
     fs_convert[23]=8
     fs_convert[24]=24
     fs_convert[25]=40
     fs_convert[26]=9
     fs_convert[27]=25
     fs_convert[28]=41
     fs_convert[29]=10
     fs_convert[30]=26
     fs_convert[31]=42
     fs_convert[32]=32
     fs_convert[33]=11
     fs_convert[34]=27
     fs_convert[35]=43
     fs_convert[36]=12
     fs_convert[37]=28
     fs_convert[38]=44
     fs_convert[39]=13
     fs_convert[40]=29
     fs_convert[41]=45
     fs_convert[42]=14
     fs_convert[43]=30
     fs_convert[44]=46
     fs_convert[45]=15
     fs_convert[46]=31
     fs_convert[47]=47
     fs_convert[48]=48
     
     strips = pd.Series(strips).map(fs_convert)

     #convert front strips to energy_scale
     if energy_scale:
         lst = np.loadtxt('clbr_coef.txt',usecols=[0,1]) # [0 .. 47][..]
         #calibration function y = a*x + b
         f_clbr = lambda st,ch: lst[st-1,0]*ch + lst[st-1,1]
         channel = np.where( (strips>0)&(strips<47),
                             f_clbr(strips,channel),channel)
         
     #back strips convert function
     def bs_convert(strip,id_):
         
         column = np.where(id_% 2 ==0,
                           strip*2 - 1 + 16*id_,
                           2*((id_-1)*8 + strip) )
         
         l = len(column)
         reverse = lambda a,b,x : (a+b)*np.ones(l) -x
         
         column=np.where( (column>= 49) & (column <= 64),
                             reverse(49,64,column), column)
         column=np.where( (column>= 33) & (column <= 48),
                             reverse(33,48,column), column)
         column=np.where( (column>= 1) & (column <= 32),
                             reverse(1,32,column), column)
         column=np.where( (column>= 65) & (column <= 128),
                             reverse(65,128,column), column)
                
         column = Ser(column)
         return column
     
     #back detectors
     back_ind1 = (data[:,0]>>8) % 16
     back_ind2 = (data[:,0]>>12)% 16
     
     back_strips1 = data[:,8] / 4096 
     back_strips1 = np.where(back_strips1,back_strips1+1,0)
     back_strips1 = bs_convert(back_strips1,back_ind1)
     
     back_strips2 = data[:,10] / 4096 
     back_strips2 = np.where(back_strips2,back_strips2+1,0)
     back_strips2 = bs_convert(back_strips1,back_ind2)
     
     back_channel1 = data[:,7] % 8192
     back_channel2 = data[:,9] % 8192
     
     #time distr
     time = data[:,3]*65536 + data[:,4]
     
     #check time gaps
     time = time - time[0]
     time = np.array(time,dtype='int64')
     time_dev = np.array(np.r_[ [0],time[1:] - time[:-1] ],dtype='int64')
     len_time = len(time)
     
     def time_correct(time_stamps,diapason):
         if len(time_stamps) > 0:
             index1 = int(time_stamps[0])
             if len(time_stamps) ==1:
                 index2 = len_time
                 time[index1:] += np.ones(len_time-index1)*diapason
             else:
                 for i in range(1,len(time_stamps)):
                     index2 = time_stamps[i]
                     time[index1:] += np.ones(len_time-index1)*diapason
                     index1 = index2
                 index2 = len(time)
                 time[index1:] += np.ones(len_time-index1)*diapason 
         return time
     
     if time_corr:
         #reset of low counter
         diapason = 65536
         time_stamps = ( time_dev<0 ).nonzero()[0]
         time = time_correct(time_stamps,diapason)
               
         #reset of high counter
         diapason = 65536*65536
         time_dev = np.array(np.r_[ [0],time[1:] - time[:-1] ],dtype='int64')
         time_stamps = ( abs(time_dev)> diapason - 1000000).nonzero()[0]
         time = time_correct(time_stamps,diapason)
     
     time_dev = np.array(np.r_[ [0],time[1:] - time[:-1] ],dtype='int64')
     time_sec = (time / 1000000)
     time_min = (time_sec / 60)
     time_hours = (time_min / 60)
     time_mks = (time % 1000000)
     time_sec %= 60
     time_min %= 60
     time_hours %=24 
     
     #TOF
     tof = data[:,6]%8192
     
     frame = DF( {'id': indexes,'strip':strips,'channel':channel,
                  'time_hours':time_hours,'time_min':time_min,'time_sec':time_sec,'time_mks':time_mks,
                  'time_dlt':time_dev,
                  'time': time,
                  'b_strip 1':back_strips1,'b_strip 2':back_strips2,
                  'b_channel1':back_channel1,'b_channel2':back_channel2,
                  'beam_marker':beam_marker,
                  'tof':tof
                  },
                 index = np.arange(data.shape[0]) )
     
     if write_file:           
         frame.to_csv(filename.split('.')[1]+'.csv')
     return frame 
     
     
     
def read_files(filenames,**argv):
    
    if type(filenames) == type(pd.DataFrame()):
        return filenames
        
    #a bit of carrying 
    #read = lambda x,y=strlength: read_file(x,y)
    #primitive parser
    def func1(s):
        if '-' in s:
            fn = s.split('-')
            if len(fn) > 2:
                raise ValueError('read_files: Wrong string format')
            numb = []
            for i in fn:
                i.replace(' ','')
                numb.append( int(i.split('.')[1]) )
            numb = range(numb[0],numb[1]+1)
            for i in range(len(numb)):
                numb[i] = 'tsn.'+str(numb[i])
            return numb
        else:
            pass
    
    if '-' in filenames:
        names = func1(filenames)
    elif ',' in filenames:
        names = filenames.split(',')
    elif len(filenames) < 8 and ('tsn.' in filenames):
        names = [filenames,]
    else:
        raise ValueError('read_files: Bad string format')
        
    frame = DF([])
    
    for i in names:
        frame = pd.concat([frame,read_file(i,**argv)])
        
    return frame
            
        

def get_front_spectrs(data,energy_scale=True,tof=False): 
    #read data 
    sample = read_files(data)
    
    #tof 
    if tof:
        sample = sample[ sample['tof']>0 ]
    else:
        sample = sample[ sample['tof']==0]
    sample = sample[ sample['id']<3 ]
    spectr = sample['channel'].groupby(sample['strip'])
    
    #choose scale
    if energy_scale:
        size = 20000
    else:
        size = 8192
    #collect histograms
    list_hist = []
    for name,group in spectr:
        list_hist.append( np.histogram(group,bins = np.arange(size))[0][1:] )
    hist = np.vstack( (list_hist[0],list_hist[1]))
    for i in xrange(2,len(list_hist)):
        hist = np.vstack( (hist,list_hist[i]))  
    #sum
    sum_spectr = list_hist[0]
    for i in range(1,len(list_hist)):
        sum_spectr += list_hist[i]
    
     #visualisation
    fig = plt.figure()
    ax1 = fig.add_axes([0.125, 0.35, 0.8, 0.6])
    #ax2 = fig.add_subplot(2,1,2,sharex = ax1)
    ax2 = fig.add_axes([0.125, 0.05, 0.64, 0.25],sharex = ax1)
    x = np.arange(1,size)
    y = np.arange(1,48)
    hist = hist[:,:-1]
    a = ax1.pcolormesh(x,y,hist)
    ax2.plot(sum_spectr,'k',linestyle='steps')
    fig.colorbar(a,ax=ax1)
    axis = [1,8192,1,47]
    ax1.axis( axis )
    ax2.axis( xmax = size )
    plt.show()
    return hist,sum_spectr
   

def pos_distr(data,tof=False):
    sample = read_files(data)
    spectr = np.histogram(sample['strip'],
                              bins=np.arange(1,48))[0]
    """
    if tof:
        spectr = np.histogram(sample['strip'].ix[sample['tof']>0],
                              bins=np.arange(1,48))[0]
    else:
        spectr = np.histogram(sample['strip'].ix[sample['tof']==0],
                              bins=np.arange(1,48))[0]"""
    back_stripes = Ser(np.where(sample['b_channel1'] > sample['b_channel2'],sample['b_strip 1'],sample['b_strip 2']))
    spectr_back = np.histogram(back_stripes,bins = np.arange(1,129))[0]
    fig,(ax1,ax2) = plt.subplots(1,2)
    ax1.plot(spectr,linestyle='steps')
    ax2.plot(spectr_back,linestyle='steps')
    plt.show()
    return spectr
   
def pos_distr2d(data):
    sample = read_files(data)
    back_stripes = Ser(np.where(sample['b_channel1'] > sample['b_channel2'],sample['b_strip 1'],sample['b_strip 2']))
    spectr = np.array( pd.crosstab(sample['strip'],back_stripes)  )
    spectr = spectr[:,1:]
    
    #visualisation
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator
    #spectr = spectr[:-1,:-1]
    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.05, 0.8, 0.95])
    levels = MaxNLocator(nbins = 10).tick_values(spectr.min(),spectr.max())
    cmap = plt.get_cmap('rainbow')#('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N,clip=True)
    #y,x = np.arange(len(spectr)),np.arange(len(spectr[0]))
    #plt.pcolormesh(x,y,spectr,cmap = cmap,norm = norm)
    a=ax.imshow(spectr,cmap = cmap,norm = norm,interpolation='nearest')
    ax.invert_yaxis()
    fig.colorbar(a,ax=ax)
    plt.show()
    
    return spectr
    
    
    
def read_times(filename):
     data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,12)[:,0:5]
     indexes = data[:,0] % 16 
     strips = data[:,2] / 4096 + 1
     strips = np.where( indexes < 3, indexes*16+strips, strips)
     time = data[:,3]*65536 + data[:,4]
     #time_block = time[1:] - time[:-1]
     frame = DF( {'id': indexes,'strip':strips,'time':time},
                 index = np.arange(data.shape[0]) )
     return frame 
     
     
     
def tof_distr(filename):
    data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,14)[:,0:7]
    print 'data_len: ',len(data)
    print '0-signals ',(data[:,1]%8192 == 0).sum()
    bit_2 = (data[:,0] >> 5) % 2
    tof = data[:,6] % 8192   
    #print tof
    tof_count = (tof > 0).sum()
    print 'tof_len ',tof_count
    print 'tof_zero_len ',(data[:,6] ==0).sum()
    print 'bit_count ', (bit_2 > 0).sum()
    tof_bit_count = ( (tof > 0)&(bit_2 == 1) ).sum()
    print 'tof&bit_count/tof_count: ',float(tof_bit_count)*100/tof_count
    bit_2_count = ( bit_2 > 0).sum()
    print 'tof&bit_count/bit_2_count: ',float(tof_bit_count)*100/bit_2_count
    x = np.arange(8192)
    distr =  np.histogram(data[:,6],bins = x)[0]
    print 'len',len(distr),len(x)
    plt.plot(x[1:-1], distr[1:], linestyle = 'steps')
    plt.show()
    
    
    
def sinc_times(filename):
     data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,12)
     indexes = data[:,0] % 16 
     #strips = data[:,2] / 4096 + 1
     #strips = np.where( indexes < 3, indexes*16+strips, strips)
     time = data[:,3]*65536 + data[:,4]
     time_sinc = data[:,11]
     #metka = (data[:,0] % 256) >> 4 
     #time_block = time[1:] - time[:-1]
     frame = DF( {'id': indexes,#'strip':strips,
                  'time':time,
                  #'metka':metka,
                  'time_sinc':time_sinc},
                 index = np.arange(data.shape[0]) )
     return frame 
     
     
     
def time_sinc_distr(filename):
    sample = sinc_times(filename)
    sample['id_dev'] = sample['id'].diff()
    sample['time_dev'] = sample['time'].diff()
    sample['sinctime_dev'] = sample['time_sinc'].diff()
    sample = sample[ (sample['id_dev'] > 0)&(sample['time_dev']<50) ]  
    maxvalue = sample['sinctime_dev'].max()
    #hist1 = np.histogram(sample['sinctime_dev'],bins = np.arange(maxvalue))
    import matplotlib.pyplot as plt
    plt.title('tsn.605 16 first + 16 third.(GTA)')
    plt.hist(sample['sinctime_dev'],np.arange(maxvalue),alpha = 0.4, facecolor = 'g')
    #plt.plot(hist1[1][1:],hist1[0],linestyle = 'steps',alpha = 0.7)
    plt.show()
"""    
def show_front(filename,strip=1):
    frame = read_files(filename)
    spectr
"""
if __name__ == '__main__':
    #print read_times('tsn.562')
    #time_sinc_distr('tsn.605')
    #frame = read_files('tsn.611-tsn.615')
    frame1 = read_files('tsn.459',strlength=14,energy_scale=True)
    #pos_distr(frame1)
    hist,sum_spectr = get_front_spectrs(frame1)
    #tof_distr('tsn.371')
    #frame2 = read_files('tsn.612,tsn.613')
    #time_sinc_distr('tsn.606')
	
