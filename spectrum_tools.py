# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 17:51:17 2014

Some handy tools based on numpy and ROOT package, wrapped into handy form.
Include:
    1. Determing background : background
    2. Smoothing.
    3. Peak searching.
@author: eastwoodknight
"""

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from read_files import get_front_spectrs



def GetHist(hist,bins):
    """get data from TH1F root hist in a form of np.array"""
    new_hist = []
    for i in xrange(bins):
        a = hist.GetBinContent(i)
        new_hist.append(a)
    return np.array(new_hist)
  
  
  
def SetHist(hist,np_hist):
    """set data from TH1F root hist in a form of np.array"""
    for i in xrange(len(np_hist)):
        hist.SetBinContent(i,np_hist[i])



def background(y,niter=30,parameters=''):
    """
    return the background of the given spectrum y;
    look into ROOT's TSpectrum manual to configure parameters"""
    lenY = len(y)
    hist = rt.TH1F('hist','Spectrum',lenY,1,lenY+1)
    SetHist(hist,y)
    
    s = rt.TSpectrum()
    sm_hist = s.Background(hist,niter,parameters)
    return GetHist(sm_hist,lenY)



def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
        
     This method is based on the convolution of a scaled window with the signal.
     The signal is prepared by introducing reflected copies of the signal 
     (with the window size) in both ends so that transient parts are minimized
     in the begining and end part of the output signal.
      
      input:
           x: the input signal 
           window_len: the dimension of the smoothing window; should be an odd integer
           window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
               flat window will produce a moving average smoothing.
   
      output:
           the smoothed signal
           
      example:
   
       t=linspace(-2,2,0.1)
       x=sin(t)+randn(len(t))*0.1
       y=smooth(x)
       
     see also: 
       
     numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
     scipy.signal.lfilter  
    """
    if x.ndim != 1:
         raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
         raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
         raise ValueError, "Too small window length."
 
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
   
   
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
         w=np.ones(window_len,'d')
    else:
         w=eval('np.'+window+'(window_len)') 
    y=np.convolve(w/w.sum(),s,mode='valid')
    y = y[(window_len/2) : -1]
    #zer = np.zeros(window_len/2)
    #y = np.r_[y,zer]
    
    return y
  
  
def search_peaks(data_x,data_y,sigma = 30,wavelet='bartlett',threshold = 0.10):
    """
    Based on scipy.signal.find_peaks_cwt peak-finding tool.
    threshold parameter allow to cut spectrum on the given level.
    Return 2 arrays, containing x and y coordinated of found peaks. 'ricker'
    """
    print 'bartlett'
    f_wavelet = eval('signal.'+wavelet)
    peak_id = signal.find_peaks_cwt(data_y,np.arange(1,sigma),wavelet=f_wavelet)
    peak_id = np.array(peak_id) - 1
    xpeaks,ypeaks=[],[]
    for i in peak_id:
        i_sample = data_y[i-3:i+3] 
        sample_max_id = 0
        for j in xrange(len(i_sample)):
            if i_sample[j] > i_sample[sample_max_id]:
                sample_max_id = j
        sample_max_id += i-3
        xpeaks.append(data_x[sample_max_id])
        ypeaks.append(data_y[sample_max_id])                     
    
     
    xpeaks,ypeaks = np.array(xpeaks),np.array(ypeaks)
    numbers = ypeaks>(ypeaks.max()*threshold)
    ypeaks = ypeaks[numbers]
    xpeaks = xpeaks[numbers]-1
    return xpeaks, ypeaks

    
"""
    MAIN TEST
"""    
if __name__ == '__main__':
    
    """
    xsample = np.arange(1,101)
    def peak(x,x0,sigm):
        return 200*np.exp( - ((x-x0)/sigm)**2/2 )
    sample = peak(xsample,30,3) + peak(xsample,47,2) + peak(xsample,84,3)+ np.random.randint(6,35,len(xsample))
    fig,ax = plt.subplots(2,1)
    ax[0].plot(sample,linestyle='steps')   
    """
    
    """ INIT"""
    ind = 30
    sample,sum_spectr = get_front_spectrs('tsn.456-tsn.458',visualize = False)
    """MAKE AXES"""
    fig,ax = plt.subplots(2,1,sharex=True)
    ax[0].plot(sample[ind],linestyle='steps') 
    """PROCESSING"""
    hists = []
    xmin, xmax = 1700, 3600
    for i in sample:
        hists.append( i[xmin:xmax] )
    sample = hists
    sample[ind] = smooth(sample[ind],9,'hanning')
    sample_background = background(sample[ind],parameters='BACK1_ORDER8,BACK1_INCLUDE_COMPTON')
    sample[ind] -= sample_background
    sample[ind][ sample[ind]<0 ] = 0
    #ax[1].plot(sample[ind],linestyle='steps',color='b') 
    hist = smooth(sample[ind],5,'bartlett')
    hist = sample[ind]
    #xsample = np.arange(xmin,xmax)
    xsample = np.arange(xmin,xmin+len(hist))
    xpeaks,ypeaks = search_peaks(xsample,hist,threshold=0.35)
    #collecting the the largest 11 peaks
    dic_peaks = dict( zip(ypeaks,xpeaks) )
    dict1 = {}
    if len(dic_peaks)>=11: 
        count_peaks = 11
    else:
        count_peaks = len(dic_peaks)
    for i in sorted(dic_peaks.keys())[-1:-count_peaks:-1]:
        dict1[i] = dic_peaks[i]
    #ypeaks,xpeaks = dict1.keys(),dict1.values()
    #sorting by xpeaks-column
    dict2 = {k: v for v,k in dict1.iteritems() }
    xpeaks,ypeaks = [],[]
    for k,v in sorted(dict2.iteritems() ):
        xpeaks.append(k)
        ypeaks.append(v)
        
    """selecting valid peaks"""
    print 'Xpeaks: ',xpeaks
    
    def get_numb(x,list1):
        for i,j in enumerate(list1):
            if j == x:
                return i
    # missing of 2nd valid from left            
    if (xpeaks[2]-xpeaks[1])> (xpeaks[3]-xpeaks[2]):
        print 'YES'
        xstart,xstop = get_numb(xpeaks[0],xsample),get_numb(xpeaks[1],xsample)
        print xstart, xstop
        if xstop - xstart != 0:
            extra_xpeaks,extra_ypeaks = search_peaks(xsample[xstart:xstop],\
                    hist[xstart:xstop],threshold = 0.12)
            print '!!',extra_xpeaks,extra_ypeaks
            if len(extra_xpeaks)>=3:
                middle_peak_num = get_numb(sorted(extra_ypeaks)[-3], extra_ypeaks)
                xpeaks.insert(1,extra_xpeaks[middle_peak_num])
                ypeaks.insert(1,extra_ypeaks[middle_peak_num])
    # deleting 2 needless peaks after 3nd
    xpeaks.pop(3)
    ypeaks.pop(3)
    xpeaks.pop(3)
    ypeaks.pop(3)
    # deleting of 4 needless peak from right
    if (len(xpeaks) >= 5) and \
        (xpeaks[-4] > xpeaks[-3] - 0.75*(xpeaks[-2]-xpeaks[-3]) ):
        xpeaks.pop(-4)
        ypeaks.pop(-4)
    """OUTPUT"""
    print xpeaks, ypeaks
    ax[1].plot(xsample,hist,linestyle='steps',color='r',linewidth=3) 
    ax[1].plot(xpeaks,ypeaks,'kD')
    plt.plot()