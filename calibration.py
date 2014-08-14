# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 15:08:07 2014

Code for processing automatic calibration of alpha-scale in our experiments.
Based on spectrum_tools.py 
 
@author: eastwoodknight
"""

from spectrum_tools import *
from read_files import get_front_spectrs
from scipy.optimize import leastsq

def calibrate_spectr(sample,xmin,xmax,visualize=True):  
    """ Return xpeaks, solution (result of linear fitting of spectrum)"""
    """MAKE AXES"""
    if visualize:
        fig,ax = plt.subplots(3,1,sharex=True)
        ax[0].set_title('Raw spectrum')
        ax[0].plot(sample,linestyle='steps') 
    
    sample = sample[xmin:xmax]    
    sample = smooth(sample,9,'hanning')
    sample_background = background(sample,parameters='BACK1_ORDER8,BACK1_INCLUDE_COMPTON')
    sample -= sample_background
    sample[ sample<0 ] = 0
    #ax[1].plot(sample[ind],linestyle='steps',color='b') 
    hist = smooth(sample,5,'bartlett')
    hist = sample
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
    def get_numb(x,list1):
        for i,j in enumerate(list1):
            if j == x:
                return i
    # deleting of 4 needless peak from right
    if (len(xpeaks) >= 5) and \
        (xpeaks[-4] > xpeaks[-3] - 0.75*(xpeaks[-2]-xpeaks[-3]) ):
        xpeaks.pop(-4)
        ypeaks.pop(-4)
    # missing of 7th valid from right
     # deleting 2 needless peaks after 3nd
    while xpeaks[-6] > xpeaks[-5] - 1.3*(xpeaks[-4]-xpeaks[-5]):
        xpeaks.pop(-6)
        ypeaks.pop(-6)
    #xpeaks.pop(-6)
    #ypeaks.pop(-6)
    if (xpeaks[-6]-xpeaks[-7])> 0.75*(xpeaks[-5]-xpeaks[-4]):
        print 'YES'
        xstart,xstop = get_numb(xpeaks[-7],xsample),get_numb(xpeaks[-6],xsample)
        print xstart, xstop
        if xstop - xstart != 0:
            extra_xpeaks,extra_ypeaks = search_peaks(xsample[xstart:xstop],\
            hist[xstart:xstop],threshold = 0.12)
            print '!!',extra_xpeaks,extra_ypeaks
            if len(extra_xpeaks)>=3:
                middle_peak_num = get_numb(sorted(extra_ypeaks)[-3], extra_ypeaks)
                xpeaks.insert(1,extra_xpeaks[middle_peak_num])
                ypeaks.insert(1,extra_ypeaks[middle_peak_num])
    #deleting all other needless peaks
    while len(xpeaks) > 8:
        xpeaks.pop(0)
        ypeaks.pop(0)
    """            
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
    """
           
    """Fitting by line"""
    if len(xpeaks) < 8:
        print 'Not enought peaks'
        print xpeaks,ypeaks
        return
        
    def residuals(coef,y,x):
        return y - coef[0]*np.ones(len(x)) - coef[1]*x
   
    p0 = (0,2) #init coefficients    
    energies = np.array([6040,6143,6264,6899,7137,7922,8699,9265])
    xpeaks = np.array(xpeaks)
    solution = leastsq(residuals,p0, args = (energies,xpeaks) )[0]
    
    """OUTPUT"""
    print 'Result (xpeaks,ypeaks):', xpeaks, ypeaks
    print 'Solution:', solution
    if visualize:
        #print len(xpeaks),len(energies)
        ax[1].set_title('Processed spectrum with marked calibration peaks')
        ax[1].plot(xsample,hist,linestyle='steps',color='b',linewidth=3) 
        ax[1].plot(xpeaks,ypeaks,'ro',linewidth=4)
        ax[2].set_title('Calibration function')
        ax[2].plot(xpeaks,energies,'ro',linewidth=4,label='Data points')
        ax[2].plot(xpeaks,solution[0]*np.ones(len(xpeaks))+solution[1]*xpeaks,linewidth=2,label='Fitting line')
        #print solution[0]*np.ones(len(xpeaks))+solution[1]*xpeaks
        ax[2].legend(loc='lower right')
        plt.show()
    return xpeaks,solution
    
if __name__ == '__main__':
    """ INIT"""
    ind = 34
    sample,sum_spectr = get_front_spectrs('tsn.456-tsn.458',energy_scale=False,visualize = False)
    
    """PROCESSING"""
    xmin, xmax = 1700, 3600
    #for i in sample:
    #   hists.append( i[xmin:xmax] )
    for ind in arange(18,27):
        hist = sample[ind]   
        calibrate_spectr(hist,xmin,xmax,visualize=True)