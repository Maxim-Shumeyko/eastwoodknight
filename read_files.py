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
from mayavi import mlab
from matplotlib import cm
import read_american_format_cy as read_am
#import calibration
#import Visualisator


def get_calibration_coefficients(filename):
    """
    Return coefficients from file: array 2xN
    coefs[0] - A, coefs[1] - B
    """
    f = open(filename,'r')
    s = f.read()[1:]
    s = s.split('\n\n')
    indexes, coefs= [],[]
    for line in s:
        line = line.split('\n')[0].split()
        indexes.append(int(line[0]))
        coefs.append( [float(line[3]),float(line[7])] )
    coefs = np.array(coefs).T #pd.DataFrame(np.array(coefs),index=indexes)
    f.close()
    return coefs
	

def read_file(filename,strlength = 14,write_file=False,energy_scale=False,
              time_corr=False,strip_convert=False, clbr_front_filename='/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_front.txt',\
			 clbr_back_filename='/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_back.txt'):

     data = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '').reshape(-1,strlength)#[:,0:9]
     beam_marker = (data[:,0]>>4)% 16
     
     #front detectors
     indexes = data[:,0] % 16 
     strips = np.array(data[:,2] / 4096 +1,dtype='uint16')
     strips = np.where( indexes < 3, indexes*16+strips, strips)
     channel = data[:,1] % 8192
     Fchannel = data[:,2]%4096
     
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
     if strip_convert or energy_scale:
         strips = np.where(indexes<3,pd.Series(strips).map(fs_convert),strips)

     #convert front strips to energy_scale
     if energy_scale:
         coefs = get_calibration_coefficients(clbr_front_filename) # 2x47
         #print 'coefs',coefs
         #calibration function y = a*x + b
         f_clbr = lambda st,ch: coefs[0][st-1]*ch + coefs[1][st-1]
         #print f_clbr(np.array([10,18,46,12]),np.array([673,505,596,601]))
         strips = np.where( strips<=47,strips,np.zeros(len(strips),dtype=np.int16))
         channel = np.where( strips&channel,f_clbr(strips,channel),np.zeros(len(strips)) )
         
     #back strips convert function     
     def bs_convert1(strip,id_):
         a = pd.read_table('/home/eastwood/codes/Python_Idle/data_processing/StripOrder.txt',sep='\s+',skiprows=53,header=None,nrows=64)
         index = a[0]
         a = pd.Series(a[1])
         a.index = index
         f = lambda x: a[x]
         strip = np.where(id_>0,strip+(id_/2)*16,np.zeros(len(id_)) )
         strip = np.array(f( strip ),dtype=np.int16)
         return strip  

     def bs_convert2(strip,id_):
         a = pd.read_table('/home/eastwood/codes/Python_Idle/data_processing/StripOrder.txt',sep='\s+',skiprows=120,header=None,nrows=64)
         index = a[0]
         #print ( id_ == 4).sum()
         a = pd.Series(a[1])
         a.index = index
         f = lambda x: a[x]
         strip = np.where(id_>0,strip+(id_/2-1)*16,np.zeros(len(id_)) )
         strip = np.array(f( strip ),dtype=np.int16)
         return strip  
         
     #reading id_s, amplitudes and strip numbers for back detectors
     back_ind1 = (data[:,0]>>8) % 16 
     back_ind2 = (data[:,0]>>12)% 16
     
     back_strips1 = data[:,8] / 4096 +1
     if strip_convert or energy_scale:
         back_strips1 = bs_convert1(back_strips1,back_ind1)
     
     back_strips2 = data[:,10] / 4096 +1 
     back_strips2 = np.where(back_strips2,back_strips2,0)
     if strip_convert or energy_scale:
         back_strips2 = bs_convert2(back_strips2,back_ind2)
     
     back_channel1 = data[:,7] % 8192
     back_channel2 = data[:,9] % 8192
     
     #convert amplitudes to energies
     if energy_scale:
#         f = open('/home/eastwood/codes/Python_Idle/data_processing/clbr_coef_back.txt','r')
#         lines = f.readlines()
#         f.close()
#         coefs = []
#         strip_list = [] # numbers of strips corresponding with existing coefs
#         for i in lines[1::11]:
#             val = i.split()
#             coefs.append( [float(val[3]),float(val[7])] )
#             strip_list.append(float(val[0]))
#         coefs = np.array(coefs).T
#         #calibration function y = a*x + b
#         f_clbr = lambda st,ch: coefs[0][st-1]*ch + coefs[1][st-1]
         coefs = get_calibration_coefficients(clbr_back_filename)
         back_strips1 = np.where( back_strips1<=128,back_strips1,np.zeros(len(back_strips1),dtype=np.int16))
         back_strips2 = np.where( back_strips2<=128,back_strips2,np.zeros(len(back_strips2),dtype=np.int16))
         #print back_strips1.dtype, back_strips2.dtype
         back_channel1 = np.where( back_channel1,f_clbr(back_strips1,back_channel1),np.zeros(len(back_channel1),dtype=np.int16) )
         back_channel2 = np.where( back_channel2,f_clbr(back_strips2,back_channel2),np.zeros(len(back_channel1),dtype=np.int16) )
                                   
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
     
     frame = pd.DataFrame( {'id': indexes,'strip':strips,'channel':channel,'Fchannel':Fchannel,
                  'time_hours':time_hours,'time_min':time_min,'time_sec':time_sec,'time_mks':time_mks,
                  'time_dlt':time_dev,
                  'time': time,
                  'synchronization_time':synchronization_time,
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
        
    frame = pd.DataFrame([])
    
    for i in names:
        frame = pd.concat([frame,read_file(i,**argv)],ignore_index=True)
        
    return frame
            
            
        
def visualize_spectrum(hist,sum_spectr):#,window=None):
    #matplotlib visualization
    
    fig = plt.figure()
    ax1 = fig.add_axes([0.125, 0.35, 0.8, 0.6])
    ax2 = fig.add_axes([0.125, 0.05, 0.64, 0.25],sharex = ax1)
    x = hist.index
    y = hist.columns
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
    """ energy_scale=True,
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
    if tof:
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
    
    visualize_spectrum(hist,sum_spectr)
            
    return hist,sum_spectr
        
        
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


def get_back_spectrs(data,tof=False,threshold=0.04,visualize=True,**argv):
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
    
    #tof 
    if tof:
        sample = sample[ sample['tof']>0 ]
    else:
        sample = sample[ sample['tof']==0]
    spectr1 = sample['b_channel1'].groupby(sample['b_strip 1']) #odd numbers of strips
    spectr2 = sample['b_channel2'].groupby(sample['b_strip 2']) #even numbers of strips
    del sample
    if len(spectr1) and len(spectr2) == 0:
        raise ValueError('No such events or empty file')
    
    #choose scale
    xsize = np.arange(8192)
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
    for (name1,group1) in spectr2:
        #collecting histograms from odd strips
        if float(name1) % 2 == 1:
            names.append(name1)
            #if 'energy_scale' in argv: #sum by window to compress spectrums
            #    if argv['energy_scale']:
            #        group1 = window_sum(group1,window )
            list_hist.append( np.histogram(group1,bins = xsize)[0][1:] )
            
    for (name2,group2) in spectr1:
        if (float(name2) % 2 == 0)&(float(name2)!=0):
            names.append(name2)
            #if 'energy_scale' in argv:
            #    if argv['energy_scale']:
            #        group2 = window_sum(group2,window )
            list_hist.append( np.histogram(group2,bins = xsize)[0][1:] )
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
    print sample
    #tof 
    if tof:
        sample = sample[ sample['tof']>0 ]
    else:
        sample = sample[ sample['tof']==0]
    spectr1 = sample['channel'].groupby(sample['b_strip 1']) #odd numbers of strips
    #spectr2 = sample['channel'].groupby(sample['b_strip 2']) #even numbers of strips
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
        print name,group
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
	set_side1 = msv_dt  & (frame['id']==3) & (frame['strip'] == strip)
	set_side1 &= np.roll(frame['id']<3,1)  
	set_front1 = np.roll(set_side1,-1)   
   
	#find synchronized pairs side-focal
	set_front2 = msv_dt & (frame['id']<3) 
	set_front2 &= np.roll(frame['id']==3,1)
	set_side2 = np.roll(set_front2,-1)
	return set_side1,set_front1, set_side2,set_front2
				

def get_side_calibration_spectrs(frame,side_coefs = get_calibration_coefficients('clbr_coef_side.txt'),front_coefs = get_calibration_coefficients('clbr_coef_front.txt'),step=10,strip = 0):
    energy_spectrs = []
    channel_spectrs = []
     
    #find pairs of front-side and side-front events
    set_side1,set_front1,set_side2,set_front2 = get_focal_side_indexes(frame,strip+11)
     
    	# read side calibration coefficients
    ff_clbr = lambda st,ch: front_coefs[0][st-1]*ch + front_coefs[1][st-1]
    fs_clbr = lambda st,ch: side_coefs[0][st]*ch + side_coefs[1][st]
#    def fs_clbr(st,ch):
#        print st-1, len(side_coefs[0])
#        return side_coefs[0][st-1]*ch + side_coefs[st-1]
     
    #group data by strip's numbers and make histogramms 
    grouped = frame[set_side1].groupby('strip')
    for i, subset in grouped:
        side_strip = int(i)
        if (side_strip>10) & (side_strip<17):
            #convert amplitudes of side detectors to energy scale	
            side_strip -= 11  
            side_energies = fs_clbr(side_strip, subset['channel'])    
            indexes = np.array(subset.index) - 1 #indexes of focal events     
            strips = np.array(frame.ix[indexes]['strip'],dtype=np.int16)
            channels = np.array(frame.ix[indexes]['channel'])
            energies = ff_clbr(strips,channels)+side_energies		     
            new_channels = (np.array(energies) - side_coefs[1][side_strip])/side_coefs[0][side_strip]    
            energy_spectrs.append(np.histogram(energies, bins = np.arange(1,20000,step)) )
            channel_spectrs.append(np.histogram(new_channels, bins = np.arange(1,20000,step)) )	
            
    grouped = frame[set_side2].groupby('strip')
    for i, subset in grouped:
        side_strip = int(i)
        if (side_strip>10) & (side_strip<17):
            #convert amplitudes of side detectors to energy scale	
            side_strip -= 11
            side_energies = fs_clbr(side_strip, subset['channel'])
            indexes = np.array(subset.index) + 1 #indexes of focal events 
            strips = np.array(frame.ix[indexes]['strip'],dtype=np.int16)
            channels = np.array(frame.ix[indexes]['channel'])
            energies = ff_clbr(strips,channels)+side_energies
            new_channels = (np.array(energies) - side_coefs[1][side_strip])/side_coefs[0][side_strip]
            energy_spectrs.append(np.histogram( energies, bins = np.arange(1,20000,step)) )
            channel_spectrs.append(np.histogram( new_channels, bins = np.arange(1,20000,step)) )
            
    #summarize histograms of f-s and s-f events by strips
    xsample = channel_spectrs[0][1][1:]
#    ch_spectrs = []
#    en_spectrs = []
#    for i in xrange(6):
#        ch_spectrs.append(np.array(channel_spectrs[0][0]) + np.array(channel_spectrs[6][0]))
#        en_spectrs.append(np.array(energy_spectrs[0][0]) + np.array(energy_spectrs[6][0]))
#        del channel_spectrs[0]
#        del energy_spectrs[0]
    channel_spectrs = np.array(channel_spectrs[0][0]) + np.array(channel_spectrs[1][0])
    energy_spectrs = np.array(energy_spectrs[0][0]) + np.array(energy_spectrs[1][0])
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
    return get_front_spectrs(data,energy_scale,tof,threshold,visualize,id_mark='<3',type_scale='Fchannel',**argv)  

def pos_distr(data,tof=False,**argv):
    sample = read_files(data, **argv )
    
    if tof:
        mask = (np.ones(len(sample),dtype=bool_))
    else:
        mask = (sample['tof']==0)
        
    hist_f = np.histogram(sample['strip'][mask],bins = np.arange(1,49))
    hist_b1 = np.histogram(sample['b_strip 1'][mask],bins = np.arange(0,129))
    hist_b2 = np.histogram(sample['b_strip 2'][mask],bins = np.arange(0,129))
    
    fig, ax = plt.subplots(2,2)
    ax[0,0].plot(hist_f[0],linestyle='steps')
    ax[0,1].plot(hist_b1[0],linestyle='steps')
    ax[1,1].plot(hist_b2[0],linestyle='steps')
        
    plt.show()
    
    return hist_f,hist_b1,hist_b2
   
def pos_distr2d(data):
    sample = read_files(data)
    back_stripes = pd.Series(np.where(sample['b_channel1'] > sample['b_channel2'],sample['b_strip 1'],sample['b_strip 2']))
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
    

def find_chaines(data,dt):
    frame = read_files(data)
    frame = frame[ 
        (frame['beam_marker']==0) &
        (frame['tof']==0) ]    
    pos_distr = frame.groupby( 
        [frame['strip'], frame['b_strip 1']] )   
    for name, group in pos_distr:
        time_group = group['time'].diff()
        group = group[ time_group < dt ]
        outputf(group,name)
        
    
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
				
				
		

    
if __name__ == '__main__':
    print 'is run'
    import os
#    os.chdir(r'./exp_data')
#    fig,ax = plt.subplots(2,1)
