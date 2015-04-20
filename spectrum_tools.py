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
         raise ValueError, "spectrum.tools smooth: Only 1-dimension arrays are valid."
    if x.size < window_len:
         raise ValueError, "spectrum.tools smooth: Input vector should be bigger than window size."
    if window_len<3:
         raise ValueError, "spectrum.tools smooth: Window length is too small ."
 
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
   
   
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
         w=np.ones(window_len,'d')
    else:
         w=eval('np.'+window+'(window_len)') 
    y=np.convolve(w/w.sum(),s,mode='valid')
    y = y[(window_len/2) :-(window_len/2)]# -1]
    #zer = np.zeros(window_len/2)
    #y = np.r_[y,zer]
    
    return y
  
  
def search_peaks(data_x,data_y,parameters): #bartlett hanning
    """
    Based on scipy.signal.find_peaks_cwt peak-finding tool.
    threshold parameter allow to cut spectrum on the given level.
    Return 2 arrays, containing x and y coordinated of found peaks. 'ricker'
    """
    data_x, data_y = np.array(data_x), np.array(data_y)
    f_wavelet = eval('signal.'+parameters.wavelet)
    peak_id = signal.find_peaks_cwt(data_y,widths=parameters.widths,wavelet=f_wavelet,min_length=parameters.min_length,\
                                    min_snr=parameters.min_snr, noise_perc=parameters.noise_perc)
    #print 'peak_id',peak_id
    if len(peak_id)==0:
        raise Exception('No peak was founded')
    peak_id = np.array(peak_id) - 1
    xpeaks,ypeaks = data_x[peak_id], data_y[peak_id]
    numbers = ypeaks>(ypeaks.max()*parameters.noise_perc)
    ypeaks = ypeaks[numbers]
    xpeaks = xpeaks[numbers]#-1
    return xpeaks, ypeaks



class Search_peak_error(Exception): pass
    
def search_peaks1(data_x,data_y,parameters): #bartlett hanning
    """
    Based on scipy.signal.find_peaks_cwt peak-finding tool.
    threshold parameter allow to cut spectrum on the given level.
    Return 2 arrays, containing x and y coordinated of found peaks. 'ricker'
    """
    data_x, data_y = np.array(data_x), np.array(data_y)
    f_wavelet = eval('signal.'+parameters.wavelet)
    peak_id = signal.find_peaks_cwt(data_y,widths=parameters.widths,wavelet=f_wavelet,min_length=parameters.min_length,\
                                    min_snr=parameters.min_snr, noise_perc=parameters.noise_perc)
    #print 'peak_id',peak_id
    if len(peak_id)==0:
        raise Search_peak_error,'No peak was founded'
    peak_id = np.array(peak_id) - 1
    xpeaks,ypeaks = data_x[peak_id], data_y[peak_id]
    numbers = ypeaks>(ypeaks.max()*parameters.noise_perc)
    ypeaks = ypeaks[numbers]
    xpeaks = xpeaks[numbers]#-1
    return xpeaks, ypeaks
#    print wavelet
#    f_wavelet = eval('signal.'+wavelet)
#    peak_id = signal.find_peaks_cwt(data_y,np.arange(1,sigma),wavelet=f_wavelet)
#    peak_id = np.array(peak_id) - 1
#    xpeaks,ypeaks=[],[]
#    for i in peak_id:
#        i_sample = data_y[i-3:i+3] 
#        sample_max_id = 0
#        for j in xrange(len(i_sample)):
#            if i_sample[j] > i_sample[sample_max_id]:
#                sample_max_id = j
#        sample_max_id += i-3
#        xpeaks.append(data_x[sample_max_id])
#        ypeaks.append(data_y[sample_max_id])                      
#    xpeaks,ypeaks = np.array(xpeaks),np.array(ypeaks)
#    numbers = ypeaks>(ypeaks.max()*threshold)
#    ypeaks = ypeaks[numbers]
#    xpeaks = xpeaks[numbers]#-1


    
   
if __name__ == '__main__':
    
    #go to the data folder
    import os
    os.chdir(r'./exp_data')
    
    xsample, ch_spectr, en_spectr = np.loadtxt('normal_spectrs.txt')
#    """MAKE AXES"""
    fig,ax = plt.subplots(2,1,sharex=True)
    ax[0].plot(xsample,ch_spectr,linestyle='steps') 
#    """PROCESSING"""
    class record: pass
    search_properties = record()
    search_properties.widths=np.arange(1,8)
    search_properties.wavelet= 'ricker'
    search_properties.min_length=3
    search_properties.min_snr=1.0
    search_properties.noise_perc=0.10
    xpeaks,ypeaks = search_peaks(xsample,ch_spectr,search_properties)
#   OUTPUT
    print xpeaks, ypeaks
    ax[1].plot(xsample,ch_spectr,linestyle='steps',color='r',linewidth=3) 
    ax[1].plot(xpeaks,ypeaks,'kD')
    plt.plot()