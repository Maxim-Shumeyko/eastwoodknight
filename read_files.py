# -*- coding: utf-8 -*-
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
    
4.  get_front_spectrs(filenames) - to get alpha front spectrs
    get_back_spectrs(filenames) 
    get_side_spectrs(filenames)
    get_front_fission_spectrs(filenames)
    ! All this functions take data from files into one big frame that may cause problems with memory
    ! (if the raw binary file is larger than 1.5GB )
    ! the way to overcome the problem is to iterate throw the group of such large files
    # also may take data from exsisted dataframe
    return list of spectrs (48 spectrs, 8192 channels each)
    - show histograms of distributions
    
5.  tof_distr(filename) - to get information about tof-distributions
    output text information
    
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
#from mayavi import mlab
from matplotlib import cm
import read_american_format_cy as read_am
import os
#import calibration
#import Visualisator

#all coefficients arrays start from 0, like [0..47] or [0..127]
#so all such array are to use in this way: coefs[0,strip_number-1] + coefs[1,strip_number-1]
def get_calibration_coefficients(filename):
    """
    Return coefficients from file: array 2xN
    coefs[0] - A, coefs[1] - B
    """
    f = open(filename,'r')
    s = f.read()[1:]
    s = s.split('\n\n')
    #indexes, coefs= [],[]
    coefs = np.empty((2,128))*np.nan
    for line in s:
        line = line.split('\n')[0].split()
        #indexes.append(int(line[0]))
        #coefs.append( [float(line[3]),float(line[7])] )
        if len(line) == 0: break
        ind = int(line[0])-1
        coefs[0][ind] = float(line[3])
        coefs[1][ind] = float(line[7])
    #coefs = np.array(coefs).T #pd.DataFrame(np.array(coefs),index=indexes)
    f.close()
    return coefs
    
    
def get_fission_calibration_coefficients(filename):
    """
    Return coefficients from file: array 2xN
    coefs[0] - A, coefs[1] - B
    """
    f = open(filename,'r')
    s = f.read()[1:]
    s = s.split('\n\n')
    #indexes, coefs= [],[]
    coefs = np.empty((3,128))*np.nan
    for line in s:
        line = line.split('\n')[0].split()
        #indexes.append(int(line[0]))
        #coefs.append( [float(line[3]),float(line[7])] )
        ind = int(line[0])-1
        coefs[0][ind] = float(line[3])
        coefs[1][ind] = float(line[7])
        coefs[2][ind] = float(line[11])
    #coefs = np.array(coefs).T #pd.DataFrame(np.array(coefs),index=indexes)
    f.close()
    return coefs


def read_file(filename,write_file=False,energy_scale=False,fission_energy_scale=False,
              time_corr=False,strip_convert=False, \
              clbr_front_filename='/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_front.txt',\
              clbr_back_filename='/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_back.txt',\
              clbr_side_filename='/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_side.txt',\
              clbr_fission_front_filename='/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_fission.txt',\
              clbr_fission_back_filename=None,\
              clbr_fission_side_filename=None):
     """
     This's the main function for getting data from experimental binary file. 
     Options:
         All following option should be [True/False]:
         write_file - to output tidy data in csv format, the name of the file would be *filename*.csv (name of the input file to read)
         energy_scale - use the option to apply alpha-scale calibration coefficients 
         fission_energy_scale - the same as previous but for fission scale
         time_corr - to apply a function to correct time of events and make them go sequentially
         strip_convert - to apply strip convertion to have strip numbers matching real physical strip numbers

         All following options are the names of files containing appropriate calibration coefficients:
         clbr_front_filename
         clbr_back_filename
         clbr_side_filename
         clbr_fission_front_filename
         clbr_fission_back_filename
         clbr_fission_side_filename
     Output: pandas.DataFrame
     """
     data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,14)#[:,0:9]     
     
     indexes = data[:,0] % 16 
     strips = np.array(data[:,2] / 4096 +1,dtype='uint16')
     #delete veto events
     not_veto = ~((indexes == 3)&(strips == 10))
     data = data[not_veto,]
     indexes = indexes[not_veto]
     strips = strips[not_veto]
     
     #front detectors
     beam_mark = (data[:,0]>>4)% 16
     strips = np.where( indexes < 3, indexes*16+strips, strips)
     channel = np.array(data[:,1] % 8192,dtype=np.int32)
     Fchannel = np.array(data[:,2]%4096,dtype=np.int32)
     
     #front strips convert func
     fs_convert = {}
     fs_convert[1]=1
     fs_convert[17]=2
     fs_convert[33]=3
     fs_convert[2]=4
     fs_convert[18]=5
     fs_convert[34]=6
     fs_convert[3]=7
     fs_convert[19]=8
     fs_convert[35]=9
     fs_convert[4]=10
     fs_convert[20]=11
     fs_convert[36]=12
     fs_convert[5]=13
     fs_convert[21]=14
     fs_convert[37]=15
     fs_convert[16]=16
     fs_convert[6]=17
     fs_convert[22]=18
     fs_convert[38]=19
     fs_convert[7]=20
     fs_convert[23]=21
     fs_convert[39]=22
     fs_convert[8]=23
     fs_convert[24]=24
     fs_convert[40]=25
     fs_convert[9]=26
     fs_convert[25]=27
     fs_convert[41]=28
     fs_convert[10]=29
     fs_convert[26]=30
     fs_convert[42]=31
     fs_convert[32]=32
     fs_convert[11]=33
     fs_convert[27]=34
     fs_convert[43]=35
     fs_convert[12]=36
     fs_convert[28]=37
     fs_convert[44]=38
     fs_convert[13]=39
     fs_convert[29]=40
     fs_convert[45]=41
     fs_convert[14]=42
     fs_convert[30]=43
     fs_convert[46]=44
     fs_convert[15]=45
     fs_convert[31]=46
     fs_convert[47]=47
     fs_convert[48]=48
     
     #we actually must apply strip convertion to have strip values matching calibration coefficients
     zeros_msv = np.array(np.zeros(len(strips)),dtype = np.int)
     if strip_convert or energy_scale or fission_energy_scale:
         strips = np.where(indexes<3,pd.Series(strips).map(fs_convert),strips)
         strips = np.where( strips<=48,strips,zeros_msv)
         strips = np.where((indexes==3)&(strips>10),strips-10,strips)
         
         a = pd.read_table('/home/eastwood/codes/Python_Idle/data_processing/StripOrder.txt',sep='\s+',skiprows=53,header=None,nrows=64)
         index = a[0]
         a = pd.Series(a[1])
         a.index = index
         f_even = lambda x: a[x]
         
     #convert front channels to alpha energy_scale
     if energy_scale and clbr_front_filename:
         coefs = get_calibration_coefficients(clbr_front_filename) # 2x47
         #calibration function y = a*x + b
         f_clbr = lambda st,ch: coefs[0][st-1]*ch + coefs[1][st-1] #alpha-scale calibration function
         channel = np.where( (indexes<3) * strips * channel,f_clbr(strips,channel),channel )
     
     #convert front channels to fission energy_scale
     if fission_energy_scale and clbr_fission_front_filename:
         coefs = get_fission_calibration_coefficients(clbr_fission_front_filename)
         F_clbr = lambda st,ch: coefs[0][st-1]*np.ones(len(ch)) + coefs[1][st-1]*np.exp(-ch/coefs[2][st-1])
         Fchannel = np.where( (indexes<3)*strips*Fchannel,F_clbr(strips,Fchannel),Fchannel )
         
     #convert side channels to alpha energy_scale
     if energy_scale and clbr_side_filename:
         coefs = get_calibration_coefficients(clbr_side_filename) # 2x47
         #calibration function y = a*x + b
         f_clbr = lambda st,ch: coefs[0][st-1]*ch + coefs[1][st-1] #alpha-scale calibration function
         #print f_clbr(np.array([10,18,46,12]),np.array([673,505,596,601]))
         channel = np.where( (indexes==3)*strips*channel,f_clbr(strips,channel),channel )
         channel = np.where( channel>0, channel, zeros_msv)
      
     #convert side channels to fission energy_scale
     if fission_energy_scale and clbr_fission_side_filename:
         coefs = get_calibration_coefficients(clbr_fission_side_filename)
         F_clbr = lambda st,ch: coefs[0][st-1]*ch + coefs[1][st-1]
         Fchannel = np.where( (indexes==3)*strips*Fchannel,F_clbr(strips,Fchannel),Fchannel )
     
     list_even = pd.read_table('/home/eastwood/codes/Python_Idle/data_processing/StripOrder.txt',sep='\s+',skiprows=53,header=None,nrows=64)
     #print list_even
     index = list_even[0]
     list_even = pd.Series(list_even[1])
     list_even.index = index
     f_even = lambda x: list_even[x]
     #back strips convert function     
     def bs_convert1(strip,id_,f=f_even):
         if strip_convert or energy_scale:          
             strip = np.where((id_>0)&(id_<=8),strip+(id_/2)*16,zeros_msv )
             strip = np.array(f( strip ),dtype=np.int16)
         else:
             #print 'yes even'
             #print 'unique ids',np.unique(id_)
             strip = np.where((id_>0)&(id_<=8),strip+(id_-1)*16,zeros_msv)
             #print 'unique',np.unique(strip)
         return strip 
         
     list_odd = pd.read_table('/home/eastwood/codes/Python_Idle/data_processing/StripOrder.txt',sep='\s+',skiprows=120,header=None,nrows=64)
     #print list_odd
     index = list_odd[0]
     #print ( id_ == 4).sum()
     list_odd = pd.Series(list_odd[1])
     list_odd.index = index
     f_odd = lambda x: list_odd[x]
     def bs_convert2(strip,id_,f=f_odd):    
         if strip_convert or energy_scale:
             strip = np.where((id_>0)&(id_<=8),strip+(id_/2-1)*16,zeros_msv )
             strip = np.array(f( strip ),dtype=np.int16)
         else:
             #print 'yes odd'
             #print 'unique ids',np.unique(id_)
             strip = np.where((id_>0)&(id_<=8),strip+(id_-1)*16,zeros_msv)
             #print 'unique',np.unique(strip)
         return strip  
         
     #reading id_s, amplitudes and strip numbers for back detectors
     back_ind1 = (data[:,0]>>8) % 16 
     back_ind2 = (data[:,0]>>12)% 16
     back_strips1 = data[:,8] / 4096 +1
     back_strips1 = bs_convert1(back_strips1,back_ind1)
     back_strips2 = data[:,10] / 4096 +1 
     back_strips2 = bs_convert2(back_strips2,back_ind2)
     back_channel1 = data[:,7] % 8192
     back_channel2 = data[:,9] % 8192
     
     Fback_channel1 = np.array(data[:,8] % 4096,dtype = np.int32)
     Fback_channel2 = np.array(data[:,10] % 4096, dtype = np.int32)     
     
     #convert back channel to alpha energy_scale
     if energy_scale and clbr_back_filename:
         coefs = get_calibration_coefficients(clbr_back_filename)
         back_strips1 = np.where( back_strips1<=128,back_strips1,zeros_msv)
         back_strips2 = np.where( back_strips2<=128,back_strips2,zeros_msv)
         #print back_strips1.dtype, back_strips2.dtype
         f_clbr = lambda st,ch: coefs[0][st-1]*ch + coefs[1][st-1]
         
         back_channel1 = np.where( back_channel1,f_clbr(back_strips1,back_channel1),back_channel1 )
         back_channel2 = np.where( back_channel2,f_clbr(back_strips2,back_channel2),back_channel2 )
    
     #convert back channel to fission energy_scale
     if fission_energy_scale and clbr_fission_back_filename:
         coefs = get_fission_calibration_coefficients(clbr_fission_back_filename) 
         F_clbr = lambda st,ch: coefs[0][st-1] + coefs[1][st-1]*np.exp(-ch/coefs[2][st-1])
         Fback_channel1 = np.where( Fback_channel1,F_clbr(back_strips1,Fback_channel1),Fback_channel1 )
         Fback_channel2 = np.where( Fback_channel2,F_clbr(back_strips2,Fback_channel2),Fback_channel2 )
     
                                  
     #time distr
     time = data[:,3]*65536 + data[:,4]
     synchronization_time = data[:,11]
     #check time gaps
     time = time - time[0]
     time = np.array(time,dtype='int64')
     time_dev = np.array(np.r_[ [0],time[1:] - time[:-1] ],dtype='int64')
     len_time = len(time)
     
     def time_correct(time_stamps,diapason):
         ind1,ind2,x = 0,0,0
         for i in time_stamps:
             ind2 = i
             if ind1*ind2:
                time[ind1:ind2] += diapason*x
             ind1 = i
             x += 1
         time[ind2:]+=diapason*x
         return time
     
     #time corrections is needed because of electric time counters sometimes reset
     if time_corr:
         #reset of low counter
         diapason = 65536
         time_stamps = ((time_dev<0)&(time_dev>-65536)).nonzero()[0]
         time = time_correct(time_stamps,diapason)
         
         #reset of high counter
         diapason = 65536**2
         time_dev = np.array(np.r_[ [0],time[1:] - time[:-1] ],dtype='int64')
         time_stamps = (time_dev<0).nonzero()[0] #!!! can be problem of adding too much time to the time of file
         print time_stamps
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
     tof = data[:,6]%4096#8192
     frame = pd.DataFrame( {'id': indexes,'strip':strips,'channel':channel,'fission_channel':Fchannel,
                  'time_hours':time_hours,'time_min':time_min,'time_sec':time_sec,'time_mks':time_mks,
                  'time_dlt':time_dev,
                  'time': time,
                  'synchronization_time':synchronization_time,
                  'back_strip_even':back_strips1,'back_strip_odd':back_strips2,
                  'back_channel_EvenStrip':back_channel1,'back_channel_OddStrip':back_channel2,
                  'back_fission_channel_EvenStrip':Fback_channel1,'back_fission_channel_OddStrip':Fback_channel2,
                  'beam_mark':beam_mark,
                  'tof':tof
                  },
                 index = np.arange(data.shape[0]) )
     
     frame = frame[frame['channel']>0]
     frame.index = np.arange(frame.shape[0])
     if write_file:           
         frame.to_csv(filename.split('.')[1]+'.csv')
         
     return frame 
     
     
     
def read_files(filenames,**argv):
    """
    Valid filename's patterns:
    tsn.xxx
    tsn.xxx, tsn.yyy, tsn.zzz, ...
    tsn.xxx - tsn.yyy
    tsn.xxx-tsn.yyy and etc.
    """
    if type(filenames) == type(pd.DataFrame()):
        return filenames
    
    #a bit of carrying 
    #read = lambda x,y=strlength: read_file(x,y)
    #primitive parser
    def func1(s):
        if '-' in s:
            fn = s.split('-')
            if len(fn) > 2:
                raise ValueError('read_files: Incorrect filenames')
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
        raise ValueError('read_files: Incorrect filenames')
        
    frame = pd.DataFrame([])
    
    for i in names:
        try:
            frame = pd.concat([frame,read_file(i,**argv)],ignore_index=True)
        except IOError,e:
            print e
            continue
        
    return frame
            
            
        
def visualize_spectrum(hist,sum_spectr):#,window=None):
    #matplotlib visualization
    
    fig = plt.figure()
    ax1 = fig.add_axes([0.125, 0.35, 0.8, 0.6])
    ax2 = fig.add_axes([0.125, 0.05, 0.64, 0.25],sharex = ax1)
    x = hist.index
    y = hist.columns
    y = np.r_[0,y]
    #in case of energy scale i use sum by window to compress spectr and make peaks more distinct
    xgrid, ygrid = np.meshgrid(x,y)
    a = ax1.pcolormesh(xgrid,ygrid,hist.as_matrix().T)
    h_min, h_max = min(hist.min()), max(hist.max())
    ticks = np.linspace(h_min,h_max,40)
    ax2.plot(x,sum_spectr,'k',linestyle='steps')
    fig.colorbar(a,ax=ax1,ticks=ticks)
    #axis = [x.min(),x.max(),y.min(),y.max()]
    #ax1.axis( axis )
    ax2.axis( xmax = max(x), ymax = max(sum_spectr) )
    plt.show()




def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    
def window_sum(a,n=4):
    return np.sum(rolling_window(np.array(a),n),axis=1)[::n]    

def get_front_spectrs(data,tof=False,threshold=0.04,visualize=True,xsize=8192,window=False,id_mark='<3',type_scale='channel',type_strip='strip',**argv): 
    """ 
    Get amplitude spectrs of front detectors from raw data.
        energy_scale - change the size of scale from 8192 to 20000 (it applies no calibrations!)
        tof - choose the type of include events by the TOF-mark
        threshold [0,1) - level of cutting of specturum (comparable with the max-height peak)
        visualize - show the distribution
        **argv - arguments to pass to read_files function
    !note: there're 48 front detectors with strip numbers 1-48; 
    Output: spectrs( pd.DataFrame [1..48]x[0..8191] , index - contains x-axis data),summary spectr
    """
    #read data 
    sample = read_files(data,**argv)
    
    #choose time-of-flight's subset  
    if tof == 'all':
        pass
    elif tof:
        sample = sample[ sample['tof']>0 ]
    else:
        sample = sample[ sample['tof']==0]

    exec( "sample = sample[ sample['id']"+id_mark+"]" ) #select events with needed id marks
    sample = sample[sample[type_strip]>0]
    spectr = sample[type_scale].groupby(sample[type_strip])
    del sample
    if len(spectr) == 0:
        raise ValueError('No such events or empty file')
    
    #choose scale
    ysize = len(spectr)
    #xsize=8192
    energy_scale = False
    #choose window	
    if 'energy_scale' in argv:
        if argv['energy_scale']:
            energy_scale = True
            window = 8
            xsize = 20000					
    if window:
		xsample=np.arange(1,xsize,window)
    else:
        xsample=np.arange(1,xsize)
    
    #collect histograms
    list_hist = []
    list_columns = []
    for name,group in spectr:
        list_hist.append( np.histogram(group,bins = xsample)[0][1:] )
        list_columns.append(name)
    hist = pd.DataFrame(np.array(list_hist).T,columns = list_columns)
    hist.index = xsample[1:-1]
    #sum
    sum_spectr = list_hist[0]
    for i in range(1,len(list_hist)):
        sum_spectr += list_hist[i]

    hist = hist.sort_index(axis=1)
    hist[ hist > max(hist.max())*threshold ] = max(hist.max())*threshold
    
    if visualize:
        visualize_spectrum(hist,sum_spectr)
            
    return hist,sum_spectr
        


def get_fission_spectrs(data,id_mark='<3',xsize=250000,window=500,**argv): 
    """ 
    Get amplitude spectrs of front detectors from raw data.
        energy_scale - change the size of scale from 8192 to 20000 (it applies no calibrations!)
        tof - choose the type of include events by the TOF-mark
        threshold [0,1) - level of cutting of specturum (comparable with the max-height peak)
        visualize - show the distribution
        **argv - arguments to pass to read_files function
    !note: there're 48 front detectors with strip numbers 1-48; 
    Output: spectrs( pd.DataFrame [1..48]x[0..8191] , index - contains x-axis data),summary spectr
    """
    if id_mark == 'side':
        id_mark = '==3'
    return get_front_spectrs(data,id_mark=id_mark,xsize=xsize,window=window,type_scale='fission_channel',**argv)
     
#apply some get_... function to get distributions from american format binaries
def get_am_distribution( data, f = lambda x: get_front_spectrs(x,tof=False,threshold=0.1,visualize=False),visualize=True ):
    Fr_obj = read_am.OpenAmFile(data,chunk_size=1000000)  
    frames = Fr_obj.get_frames()
    hist,sum1 = f(frames.next())
    for frame in frames:
         hist_,sum1_ = f(frame)   
         hist += hist_
         sum1 += sum1_
    if visualize:
        visualize_spectrum(hist,sum1) 
    return hist,sum1
    
    
    
def get_side_spectrs(data,energy_scale=True,tof=False,threshold=0.04,visualize=True,**argv): 
    """ 
    Get amplitude spectrs of side detectors from raw data.
        energy_scale - change the size of scale from 8192 to 20000 (it applies no calibrations!)
        tof - choose the type of include events by the TOF-mark
        threshold [0,1) - level of cutting of specturum (comparable with the max-height peak)
        visualize - show the distribution
        **argv - arguments to pass to read_files function
    !note: there're 6 side detectors with strip numbers 11-16; strip №10 - is a veto detector
    Output: spectrs( 2Darray [0..47]x[0..8191] ),summary spectr
    """
    #(data,tof=False,threshold=0.04,visualize=True,id_mark='<3',type_scale='channel',type_strip='strip',**argv): 
    return get_front_spectrs(data,tof,threshold,visualize,id_mark='==3',**argv)
"""
def get_am_side_spectrs( data,tof=False,threshold=0.04,visualize=True,**argv):
    Fr_obj = read_am.OpenAmFile(data,chunk_size=1000000)  
    frames = Fr_obj.get_frames()
    hist,sum1 = get_side_spectrs(frames.next(),tof=tof,threshold=threshold,visualize=False)
    for frame in frames:
        hist_,sum1_ = get_side_spectrs(frame,tof=tof,threshold=threshold,visualize=False)    
        hist += hist_
        sum1 += sum1_
    if visualize:
        visualize_spectrum(hist,sum1) 
    return hist,sum1
"""


def get_back_spectrs(data,tof=False,threshold=0.04,visualize=True,fission_scale=False,**argv):
    """ energy_scale=True,
    Get amplitude spectrs of back detectors from raw data.
        energy_scale - change the size of scale from 8192 to 20000 (it applies no calibrations!)
        tof - choose the type of include events by the TOF-mark
        threshold [0,1) - level of cutting of specturum (comparable with the max-height peak)
        visualize - show the distribution
        fission_scale - to use fission channels instead of alpha
        **argv - arguments to pass to read_files function
    !note: there're 128 back detectors, which separated to odd and even numbers (signals are splitted to different electric chains); 
    Output: spectrs( pd.DataFrame [1..48]x[0..8191] , index - contains x-axis data),summary spectr
    """
    #read data 
    sample = read_files(data,**argv)
    
    #tof 
    if tof:
        sample = sample[ sample['tof']>0 ]
    else:
        sample = sample[ sample['tof']==0]
    
    if fission_scale:
        channel_even = 'back_fission_channel_EvenStrip'
        channel_odd = 'back_fission_channel_OddStrip'
    else:
        channel_even = 'back_channel_EvenStrip'
        channel_odd = 'back_channel_OddStrip'
    spectr1 = sample[channel_even].groupby(sample['back_strip_even']) #odd numbers of strips
    spectr2 = sample[channel_odd].groupby(sample['back_strip_odd']) #even numbers of strips
    del sample
    if len(spectr1) and len(spectr2) == 0:
        raise ValueError('No such events or empty file')
    
    #choose scale
    xsample = np.arange(8192)
    energy_scale = False
    if 'energy_scale' in argv:
        if argv['energy_scale']:
            energy_scale = True
            window = 8
            xsize = 20000					
            xsample=np.arange(1,xsize,window)
    if fission_scale:
        xsample = np.arange(1,100000,100)
#    if 'energy_scale' in argv:
#        if argv['energy_scale']:
#            energy_scale = True
#            xsize = np.arange(20000)
#            window = 8
#            xsize = xsize[::window ]
            
    #collect histograms
    list_hist = []
#    names = []
#    for (name1,group1) in spectr2:
#        #collecting histograms from odd strips
#        if float(name1) % 2 == 1:
#            names.append(name1)
#            #if 'energy_scale' in argv: #sum by window to compress spectrums
#            #    if argv['energy_scale']:
#            #        group1 = window_sum(group1,window )
#            list_hist.append( np.histogram(group1,bins = xsample)[0][1:] )
#            
#    for (name2,group2) in spectr1:
#        if (float(name2) % 2 == 0)&(float(name2)!=0):
#            names.append(name2)
#            #if 'energy_scale' in argv:
#            #    if argv['energy_scale']:
#            #        group2 = window_sum(group2,window )
#            list_hist.append( np.histogram(group2,bins = xsample)[0][1:] )
#    hist = np.vstack( (list_hist[0],list_hist[1]))
#    for i in xrange(2,len(list_hist)):
#        hist = np.vstack( (hist,list_hist[i]))  
#    ysize = hist.shape[0]
    names = set(spectr1.groups.keys()).union(spectr2.groups.keys())
    if len(names) == 0:
        raise Exception('Empty strip"s set')
    
    new_names = []    
    for name in names:
        if float(name) <= 0 or float(name) >128:
            continue
        elif (name in spectr1.groups.keys()) and (name in spectr2.groups.keys()):
            list_hist.append( np.histogram(spectr1.get_group(name),bins = xsample)[0][1:] + np.histogram(spectr2.get_group(name),bins = xsample)[0][1:] )
            new_names.append(name)
        elif (name in spectr1.groups.keys()):
            list_hist.append( np.histogram(spectr1.get_group(name),bins = xsample)[0][1:]) 
            new_names.append(name)
        elif (name in spectr2.groups.keys()):
            list_hist.append( np.histogram(spectr2.get_group(name),bins = xsample)[0][1:]) 
            new_names.append(name)
    
    names = new_names
    hist = np.zeros((len(names),len(list_hist[0])))
    for i in xrange(len(names)):
        hist[i] = np.array(list_hist[i])
    #calculate the sum of spectrums
    sum_spectr = list_hist[0]
    for i in xrange(1,len(list_hist)):
        sum_spectr += list_hist[i]   
        
    hist[ hist > hist.max()*threshold ] = hist.max()*threshold
    hist = pd.DataFrame(hist.T,columns=names)
    hist.index = xsample[1:-1]
    hist = hist.sort_index(axis=1)
    #hist[ hist > max(hist.max())*threshold ] = max(hist.max())*threshold
    
    #visualisation
    if visualize:
        if energy_scale:
            visualize_spectrum(hist,sum_spectr)#,window=window)
        else:
            visualize_spectrum(hist,sum_spectr)
    return hist,sum_spectr
 

def get_am_back_spectrs(data,tof=False,threshold=0.04,visualize=True,**argv):
    """ energy_scale=True,
    Get amplitude spectrs of back detectors from raw data.
        energy_scale - change the size of scale from 8192 to 20000 (it applies no calibrations!)
        tof - choose the type of include events by the TOF-mark
        threshold [0,1) - level of cutting of specturum (comparable with the max-height peak)
        visualize - show the distribution
        **argv - arguments to pass to read_files function
    !note: there're 128 back detectors, which separated to odd and even numbers (signals are splitted to different electric chains); 
    Output: spectrs( 2Darray [0..47]x[0..8191] ),summary spectr
    """
    #read data 
    sample = read_files(data,**argv)
    sample = sample[sample['id']==4]
    #print sample
    #tof 
    if tof:
        sample = sample[ sample['tof']>0 ]
    else:
        sample = sample[ sample['tof']==0]
    spectr1 = sample['channel'].groupby(sample['back_strip_even']) #odd numbers of strips
    #spectr2 = sample['channel'].groupby(sample['back_strip_odd']) #even numbers of strips
    del sample
    if len(spectr1) == 0:
        raise ValueError('No such events or empty file')
    
    #choose scale
    xsize = np.arange(16384)
    energy_scale = False
    if 'energy_scale' in argv:
        if argv['energy_scale']:
            energy_scale = True
            xsize = np.arange(20000)
            window = 8
            xsize = xsize[::window ]
            
    #collect histograms
    list_hist = []
    names = []
#    for (name1,group1) in spectr2:
#        #collecting histograms from odd strips
#        if float(name1) % 2 == 1:
#            names.append(name1)
#            #if 'energy_scale' in argv: #sum by window to compress spectrums
#            #    if argv['energy_scale']:
#            #        group1 = window_sum(group1,window )
#            list_hist.append( np.histogram(group1,bins = xsize)[0][1:] )
            
    for (name,group) in spectr1:
        #print name,group
        #if (float(name2) % 2 == 0)&(float(name2)!=0):
        names.append(name)
        list_hist.append( np.histogram(group,bins = xsize)[0][1:] )
    hist = np.vstack( (list_hist[0],list_hist[1]))
    for i in xrange(2,len(list_hist)):
        hist = np.vstack( (hist,list_hist[i]))  
    ysize = hist.shape[0]
    
    #calculate the sum of spectrums
    sum_spectr = list_hist[0]
    for i in xrange(1,len(list_hist)):
        sum_spectr += list_hist[i]     
    hist = pd.DataFrame(hist.T,columns=names)
    hist = hist.sort_index(axis=1)
    hist[ hist > max(hist.max())*threshold ] = max(hist.max())*threshold
    
    #visualisation
    if visualize:
        if energy_scale:
            visualize_spectrum(hist,sum_spectr,window=window)
        else:
            visualize_spectrum(hist,sum_spectr)
    return hist,sum_spectr
    
    

def get_focal_side_indexes(frame,strip):      
	msv_dt = frame['synchronization_time'].diff().abs() < 20  
	#find synchronized pairs focale-side
	set_side1 = (msv_dt  & (frame['id']==3) & (frame['strip'] == strip)).as_matrix()
	set_side1 &= np.roll(frame['id']<3,1)  
	set_front1 = np.roll(set_side1,-1)   
   
	#find synchronized pairs side-focal
	set_front2 = (msv_dt & (frame['id']<3) ).as_matrix()
	set_front2 &= np.roll(frame['id']==3,1) & np.roll(frame['strip'] == strip,1)
	set_side2 = np.roll(set_front2,-1)  
	return set_side1,set_front1, set_side2,set_front2
				

def get_side_calibration_spectrs(frame,side_coefs = get_calibration_coefficients('clbr_coef_side.txt'),front_coefs = get_calibration_coefficients('clbr_coef_front.txt'),fission_scale=False,step=10,strip = 1):
    """
    Get focal-side alph-scale spectr for performing calibration of side detectors.
    Input:
        frame - initial dataset in pd.DataFrame format
        side_coefs,front_coefs - initial coefficients for making a figure
        step - width of bins in resulting figure, it is to make peaks more narrow
        strip - number of side strip for making a spectrs from it
        
    Output: xsample (numpy.ndarray,x-scale for figures),channel_spectrs (numpy.ndarray),energy_spectrs (numpy.ndarray)
    """
    energy_spectrs = []
    channel_spectrs = []
     
    #find pairs of front-side and side-front events
    #print strip+11
    set_side1,set_front1,set_side2,set_front2 = get_focal_side_indexes(frame,strip) # strip+11
    #print sum(set_side1),sum(set_side2)
     
    # read side calibration coefficients
    ff_clbr = lambda st,ch: front_coefs[0][st-1]*ch + front_coefs[1][st-1]
    fs_clbr = lambda st,ch: side_coefs[0][st-1]*ch + side_coefs[1][st-1]
    if fission_scale:
        ff_clbr = lambda st,ch: front_coefs[0][st-1] + front_coefs[1][st-1]*np.exp(-ch/front_coefs[2][st-1])
        #fs_clbr = lambda st,ch: side_coefs[0][st-1] + side_coefs[1][st-1]*np.exp(-ch/side_coefs[2][st-1])
        channel = 'fission_channel'
        energy_msv_size = 250000
        step_en = 200
    else:
        channel = 'channel'
        energy_msv_size = 20000
        step_en = step
#    def fs_clbr(st,ch):
#        print st-1, len(side_coefs[0])
#        return side_coefs[0][st-1]*ch + side_coefs[st-1]
     
    #group data by strip's numbers and make histogramms 
    grouped = frame[set_side1].groupby('strip')
    for i, subset in grouped:
        side_strip = int(i)
#        if (side_strip>10) & (side_strip<17):
#            #convert amplitudes of side detectors to energy scale	
#            side_strip -= 11  
        side_energies = fs_clbr(side_strip, subset[channel])    
        indexes = np.array(subset.index) - 1 #indexes of focal events     
        strips = np.array(frame.ix[indexes]['strip'],dtype=np.int16)
        channels = np.array(frame.ix[indexes][channel])
        energies = ff_clbr(strips,channels)+side_energies
        #print 'energies: ',ff_clbr(strips,channels),side_energies
        new_channels = (np.array(energies) - side_coefs[1][side_strip-1])/side_coefs[0][side_strip-1]    
        energy_spectrs.append(np.histogram(energies, bins = np.arange(1,energy_msv_size,step_en)) )
        #print '1 sum of energy msv nonzero items: ',max(side_energies),max(energies)
        channel_spectrs.append(np.histogram(new_channels, bins = np.arange(1,20000,step)) )	
            
    grouped = frame[set_side2].groupby('strip')
    for i, subset in grouped:
        side_strip = int(i)
#        if (side_strip>10) & (side_strip<17):
#            #convert amplitudes of side detectors to energy scale	
#            side_strip -= 11
        side_energies = fs_clbr(side_strip, subset[channel])
        indexes = np.array(subset.index) + 1 #indexes of focal events 
        strips = np.array(frame.ix[indexes]['strip'],dtype=np.int16)
        channels = np.array(frame.ix[indexes][channel])
        energies = ff_clbr(strips,channels)+side_energies
        #print 'energies: ',ff_clbr(strips,channels),side_energies
        new_channels = (np.array(energies) - side_coefs[1][side_strip-1])/side_coefs[0][side_strip-1]
        energy_spectrs.append(np.histogram( energies, bins = np.arange(1,energy_msv_size,step_en)) )
        #print '2 sum of energy msv items: ',max(side_energies),max(energies)
        channel_spectrs.append(np.histogram( new_channels, bins = np.arange(1,20000,step)) )
        
    #summarize histograms of f-s and s-f events by strips
    xsample = channel_spectrs[0][1][1:]
    xsample_en = energy_spectrs[0][1][1:]
#    ch_spectrs = []
#    en_spectrs = []
#    for i in xrange(6):
#        ch_spectrs.append(np.array(channel_spectrs[0][0]) + np.array(channel_spectrs[6][0]))
#        en_spectrs.append(np.array(energy_spectrs[0][0]) + np.array(energy_spectrs[6][0]))
#        del channel_spectrs[0]
#        del energy_spectrs[0]
    channel_spectrs = np.array(channel_spectrs[0][0]) + np.array(channel_spectrs[1][0])
    energy_spectrs = np.array(energy_spectrs[0][0]) + np.array(energy_spectrs[1][0])
    if fission_scale:
        return xsample,channel_spectrs,xsample_en,energy_spectrs
    return xsample,channel_spectrs,energy_spectrs#ch_spectrs,en_spectrs     

		
    
def get_front_fission_spectrs(data,energy_scale=True,tof=False,threshold=0.04,visualize=True,**argv): 
    """ 
    Get amplitude spectrs of front detectors in high-energy scale from raw data.
        energy_scale - change the size of scale from 8192 to 20000 (it applies no calibrations!)
        tof - choose the type of include events by the TOF-mark
        threshold [0,1) - level of cutting of specturum (comparable with the max-height peak)
        visualize - show the distribution
        **argv - arguments to pass to read_files function
    Output: spectrs( 2Darray [0..47]x[0..8191] ),summary spectr
    """
    return get_front_spectrs(data,energy_scale,tof,threshold,visualize,id_mark='<3',type_scale='fission_channel',**argv)  

def pos_distr(data,tof=False,**argv):
    sample = read_files(data, **argv )
    
    if tof:
        mask = (np.ones(len(sample),dtype=bool_))
    else:
        mask = (sample['tof']==0)
        
    hist_f = np.histogram(sample['strip'][mask],bins = np.arange(1,49))
    hist_b1 = np.histogram(sample['back_strip_even'][mask],bins = np.arange(0,129))
    hist_b2 = np.histogram(sample['back_strip_odd'][mask],bins = np.arange(0,129))
    
    fig, ax = plt.subplots(2,2)
    ax[0,0].plot(hist_f[0],linestyle='steps')
    ax[0,1].plot(hist_b1[0],linestyle='steps')
    ax[1,1].plot(hist_b2[0],linestyle='steps')
        
    plt.show()
    
    return hist_f,hist_b1,hist_b2
   
def pos_distr2d(data):
    sample = read_files(data)
    back_stripes = pd.Series(np.where(sample['back_channel_EvenStrip'] > sample['back_channel_OddStrip'],sample['back_strip_even'],sample['back_strip_odd']))
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

def outputf(data,name):
    data.to_csv('chain'+str(name))
    

#def find_chaines(data,dt):
#    frame = read_files(data)
#    frame = frame[ 
#        (frame['beam_mark']==0) &
#        (frame['tof']==0) ]    
#    pos_distr = frame.groupby( 
#        [frame['strip'], frame['back_strip_even']] )   
#    for name, group in pos_distr:
#        time_group = group['time'].diff()
#        group = group[ time_group < dt ]
#        outputf(group,name)
        
        
def get_total_energy(frame,scale = 'alpha'):
    """
    Функция суммирует энергии парных событий типа focal-side и записывает эту общую энергию в фокальное событие,
    содержащее координаты - номера фронтального и заднего стрипов.
    Output: frame (pd.DataFrame)
    """
    if scale == 'alpha':
        channel = 'channel'
        total_energy = 'total_energy_A'
    elif scale == 'fission':
        channel = 'fission_channel'
        total_energy = 'total_energy_F'
    set_energies = np.array(frame[channel])
    for i in np.arange(1,2): #complimented by the numbers of side detectors 
    
        set_side1,set_front1, set_side2,set_front2 = get_focal_side_indexes(frame,i) # i - number of a side detector
        index0 = np.nonzero(set_front1)[0]
        set0 = set_energies[set_front1] + set_energies[set_side1]
        set_energies[index0] = set0
        
        index1 = np.nonzero(set_front2)[0]
        set1 = set_energies[set_front2] + set_energies[set_side2]
        set_energies[index1] = set1
        
    frame[total_energy] = pd.Series(set_energies)
         
#        index0 = frame.ix[set_front1].index
#        set0 = pd.Series(frame.ix[set_front1][channel].as_matrix() + frame.ix[set_side1][channel].as_matrix(),index=index0)
#        index1 = frame.ix[set_front2].index
#        set1 = pd.Series(frame.ix[set_front2][channel].as_matrix() + frame.ix[set_side2][channel].as_matrix(),index=index1)

#        frame[total_energy] = np.where()
#        frame[total_energy] = set1
        #frame.ix[set_front1][total_energy] = frame.ix[set_front1][channel] + frame.ix[set_side1][channel]
        #frame.ix[set_front1][total_energy] = frame.ix[set_front1][channel] + frame.ix[set_side1][channel]
    
    return frame


def find_chaines(frame,time_dlt=25,tof_mark=False,energy_level_min=8000,energy_level_max=12000,output = 'chains.txt',tof_condition=True):
    """
    The function searches chains of alpha decays with given parameters.
    Input:
        frame - pd.DataFrame as a result of read_file(s) function
        time_dlt - time window to select events consequence into a chain
        energy_level_min,energy_level_max - interval of alpha-energies to find
        chain_length - maximum of chain's length
    Output:
        output.txt - file with a report 
    """
    frame = get_total_energy(frame)
    frame = frame[frame['id']!=3]
    report = ''
    for name,group in frame.groupby(['strip','back_strip_even']):
        chains = find_chaines_in_mesh(group,time_dlt=time_dlt,tof_mark=tof_mark,energy_level_min=energy_level_min,energy_level_max=energy_level_max,tof_condition=tof_condition)
        if len(chains):
            report += '\n\n\n'+str(name)
            for chain in chains:
                report += '\n' + chain.to_string()
                
    for name,group in frame.groupby(['strip','back_strip_odd']):
        chains = find_chaines_in_mesh(group,time_dlt=time_dlt,tof_mark=tof_mark,energy_level_min=energy_level_min,energy_level_max=energy_level_max,tof_condition=tof_condition)
        if len(chains):
            report += '\n\n\n'+str(name)
            for chain in chains:
                report += '\n' + chain.to_string()
    if output:            
        f = open(output,'w')
        f.write(report)
        f.close()
    else:
        return report
                
            
def find_chaines_in_mesh(frame,time_dlt=35,tof_mark=False,energy_level_min=8000,energy_level_max=12000,tof_condition=True):
    #selecting pairs of valid events and putting them into lists
    #energy_conditions = (frame['total_energy_A']>energy_level_min)& (frame['total_energy_A']<energy_level_max)
    resume = (frame['time'].diff().abs() <time_dlt) #& (frame['id']!=3)#&\
#        energy_conditions & np.roll(energy_conditions,1) 
#        ( (frame['back_strip_even']>0)  ) &\
# | (frame['back_strip_odd']>0)\
#(frame['tof']== tof_mark) &\
#( frame['time'].diff()<300) &\
#    frame = frame[resume]
#    energy_conditions = (frame['total_energy_A']>energy_level_min)& (frame['total_energy_A']<energy_level_max)
#    resume = energy_conditions & np.roll(energy_conditions,1)
    if not resume.any():
        return []
    resume = np.array(resume, dtype = np.bool)
    resume_back = np.roll(resume,-1)
    resume = np.array(frame[resume].index, dtype = np.int)
    resume_back = np.array(frame[resume_back].index, dtype = np.int)
    
    #processing - collecting chains and putting them into one list
    chains = []
    i_pre = resume[0]
    j_pre = resume_back[0]
    columns = ['synchronization_time','total_energy_A','id','strip','back_strip_even','back_strip_odd','tof','time','channel']
    chain = pd.DataFrame(frame.ix[[j_pre,i_pre]] [columns])
    if len(resume) > 1:
        for i,j in zip(resume[1:],resume_back[1:]):
            if j != i_pre:
                chain = chain[(chain['total_energy_A']>energy_level_min)& (chain['total_energy_A']<energy_level_max)]
                if len(chain)>1:
                    chains.append(chain)
                else:
                    chain = []
                chain = pd.DataFrame(frame.ix[[j,i]] [columns])
            elif j == i_pre:
                chain = pd.concat([chain,frame.ix[[i]] [columns] ])
            i_pre = i
    chain = chain[(chain['total_energy_A']>energy_level_min)& (chain['total_energy_A']<energy_level_max)]
    if len(chain)>1:
        chains.append(chain)
    
    #processing - selecting chains which сonform energy conditions
#    chains_new = []
#    for i,chain in enumerate(chains):
#        chain = chain[(chain['total_energy_A']>energy_level_min)& (chain['total_energy_A']<energy_level_max)]
#        if chain.shape[0]:
#            chains_new.append(chain)
#    
#    chains = chains_new
    
    #processing - removing all chains without tof-marked events    
    if tof_condition:
        indexes_ = []
        for i,chain in enumerate(chains):
            if not chain['tof'].any():
                indexes_.append(i)
                
        for i in reversed(indexes_):
            del chains[i]
            
    #postprocessing - making fine view of chain's frames
        
    return chains
                
    
def read_times(filename):
     data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,12)[:,0:5]
     indexes = data[:,0] % 16 
     strips = data[:,2] / 4096 + 1
     strips = np.where( indexes < 3, indexes*16+strips, strips)
     time = data[:,3]*65536 + data[:,4]
     #time_block = time[1:] - time[:-1]
     frame = pd.DataFrame( {'id': indexes,'strip':strips,'time':time},
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
     frame = pd.DataFrame( {'id': indexes,#'strip':strips,
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


def diff_time_distr(filename,xmax=100):
    frame = read_file(filename)['time']
    frame = frame.diff()
    frame = frame[frame<200]
    hist,bins = np.histogram(frame,bins=np.arange(0,xmax,5))
    plt.grid(True)
    plt.xlabel('mks')
    plt.ylabel('counts')
    plt.plot(bins[1:],hist,linestyle='steps',linewidth=3,color='b')
    plt.title('Time differences in interval 0-%d mkseconds' % xmax)
    plt.show()    

# distibution of average counts per second
#file = [...some filenames]
#k = [234.,235325.,...]	
def count_rate_front(files,times,filename):
    f = open(filename,'a')	
    for file_,coef in zip(files,times):
        frame = read_file(file_,strip_convert=True)
        frame = frame[frame['id']<3]
        a = frame.groupby('strip').count()['time']
        frame1 = pd.DataFrame({'Counts':a,'Average counts per second':a/coef})
        f.write('\n'+file_+'\n'+frame1.to_string())
        f.write('\n%24s%s%4s%s\n' %(' ',str(len(frame)/coef),' ',str(len(frame))))
        
    f.close()
    
def count_rate_back(files,times,filename):
    f = open(filename,'a')	
    for file_,coef in zip(files,times):
        frame = read_file(file_,strip_convert=True)
        frame = frame[frame['id']<3]
        a = frame.groupby('back_strip_odd').count()['time']
        frame1 = pd.DataFrame({'Counts':a,'Average counts per second':a/coef})
        f.write('\n'+file_+'\n'+frame1.to_string())
        a = frame.groupby('back_strip_even').count()['time']
        frame1 = pd.DataFrame({'Counts':a,'Average counts per second':a/coef})
        f.write('\n'+frame1.to_string())
        f.write('\n%24s%s%4s%s\n' %(' ',str(len(frame)/coef),' ',str(len(frame))))
        
    f.close()
    
			
		

    
if __name__ == '__main__':
    print 'is run'
    import os
#    os.chdir(r'./Files88_91')
#    frame = read_files('tsn.88-tsn.91',strip_convert=True)
#    os.chdir(r'./exp_data1')
#    frame = read_files('tsn.35-tsn.48',strip_convert=True)
#    hist,hsum  = get_front_spectrs(frame,tof=False,type_scale='fission_channel')
#    hist1,hsum1 = get_front_spectrs(frame,tof=True,type_scale='fission_channel')
#    
#    fig,ax = plt.subplots(2,1)
