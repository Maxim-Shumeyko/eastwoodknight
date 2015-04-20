# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 20:45:04 2015

@author: eastwood
"""

from read_files import *
from read_american_format_cy import *
from calibration import *
import numpy as np
 
def visualize_side_spectr(xsample,spectrs):
    sample1 = smooth(spectrs,7,'hanning')
    sample1[sample1<0] = 0
    sample = sample1 - background(sample1,niter=3,parameters='BACK1_ORDER8,BACK1_INCLUDE_COMPTON')
    sample[sample<0] = 0
    fig,ax = plt.subplots(3,1,sharex=True,sharey=True)		
    ax[0].plot(xsample,sample,linestyle='steps')
    ax[1].plot(xsample,sample1,linestyle='steps')
    ax[2].plot(xsample,spectrs,linestyle='steps')
    plt.show()
    #return sample
 
class Bad_result(Exception): pass
       
    
def calibrate_side_calibration_spectrs(frame,side_coefs,front_coefs,calibration_properties,calibration_properties1,filter_properties,search_properties,ind,flag,step=10):
    xsample,ch_spectrs, en_spectrs = get_side_calibration_spectrs(frame, side_coefs,front_coefs,step=step,strip=ind)
    #np.savetxt('noise_spectrs.txt',(xsample,ch_spectrs[ind],en_spectrs[ind]))        
    xpeaks,coef=calibrate_area(xsample,en_spectrs,4000,11500,calibration_properties,filter_properties,search_properties)
    #xpeaks1 = xpeaks
    #print 'xpeaks,coef', xpeaks,coef
    #print 'energy scale',make_report(xpeaks,coef,ind,energies=calibration_properties.energies)
    a = side_coefs[0,ind]*coef[1] #A coef[1]
    b = side_coefs[1,ind]*coef[1] + coef[0] #B coef[0] 
    side_coefs[0,ind] = a
    side_coefs[1,ind] = b
    control_S2 = sum((xpeaks - np.array(calibration_properties.energies))**2)**0.5/len(xpeaks)
    print ' Control: ', control_S2
    print
    if control_S2 > 350:
        raise Bad_result, 'wrong coefficients or bad calibration_properties'
    #calculate calibration boundaries of spectrums
    xmin = int((4500 - b)/a) # a == side_coefs[0][i] 
    xmax = int((10500 - b)/a)  # b == side_coefs[1][i] 
    xsample,ch_spectrs, en_spectrs = get_side_calibration_spectrs(frame, side_coefs,front_coefs,step=step,strip=ind)
    xpeaks,coef=calibrate_area(xsample,ch_spectrs,xmin,xmax,calibration_properties1,filter_properties,search_properties) 
    return coef,xpeaks, control_S2
 
   
#CONSTRUCT AND SHOW FOCAL-SIDE ENERGY SPECTRS 
    
def make_focal_side_calibration(frame,front_coefs,side_coefs,calibration_properties,calibration_properties1,search_properties,filter_properties,ind = 1,iterations=3):   
    #construct calibration objects of properties
    
    control, coefs, peaks = [],[],[]    
    side_coefs_start = np.array(side_coefs)
    coef = [side_coefs[0,ind],side_coefs[1,ind] ]
    xpeaks = False
    for i in xrange(iterations):
        try:
            coefs.append(coef)
            peaks.append(xpeaks)
            coef,xpeaks,control_S2 = calibrate_side_calibration_spectrs(frame,side_coefs,front_coefs,calibration_properties,calibration_properties1,filter_properties,search_properties,ind,flag=True,step=10)
            control.append(control_S2)
            side_coefs[0,ind] = coef[1] #A coef[1]
            side_coefs[1,ind] = coef[0] #B coef[0] 
        except (Calibration_Error,Bad_result),e:
            while len(coefs) > len(control):
                coefs.pop()
            side_coefs[0,ind] += 0.05*np.random.randn()
            side_coefs[1,ind] += 20*np.random.randn()
        print 'coefs:',side_coefs[0,ind],side_coefs[1,ind]
    if len(coefs) == 0:
        print 'Bad initial coefficients'
        return
    ind1 = np.argmin(np.array(control))   
    side_coefs[0,ind] = coefs[ind1][1] #A coef[1]
    side_coefs[1,ind] = coefs[ind1][0] #B coef[0]
    print 'Control: ', control[ind1]
    return (coefs[ind1][1], coefs[ind1][0]), peaks[ind1]
    
if __name__ == '__main__':

# READ DATA
#go to the data folder and read kamak data
    import os
    side_coefs = get_calibration_coefficients('clbr_coef_side.txt')
    #print 'right coefs:',side_coefs,'\n'
    front_coefs = get_calibration_coefficients('clbr_coef_front.txt')
    os.chdir(r'./exp_data1')
    frame = read_files('tsn.35-tsn.61',energy_scale = False, strip_convert = True) #-tsn.61

#read am file and plot front and back distributions
#    frame1 = read_am.OpenAmFile('11Feb_15_18.bin',chunk_size=10000000).get_frames().next()
#    hist,hsum = get_front_spectrs(frame1,xsize=16384,threshold=0.02)
#    hist_b,hsum_b = get_am_back_spectrs(frame1,threshold=0.02)
    
#use strip_convert to make strip numbers coincide with true physical sequence of cable connections
#use energy_scale to apply calibration coefficients (files with coefficients can be also included)	
#    frame = read_file('tsn.458',strip_convert=True,energy_scale=True)

# GET DISTRIBUTIONS	
#get distribution of back spectrs
#    hist,sum_spectr = get_back_spectrs('tsn.164',strip_convert=True,threshold=1.0,energy_scale=True)
#    x,y = hist.shape
#    xsample=np.arange(0,x)*8+1
#    ax[0].plot(xsample,hist[2],linestyle='steps')
#    ax[0].set_xlim(1,x*8+1)
#and front spectrs
#    hist,sum_spectr = get_front_spectrs(frame)#,energy_scale=True) #'energy_scale' provides energy spectr    
#    ax[1].plot(hist.index,sum_spectr,linestyle='steps')
#    plt.show()
				
##get and show focal-side energy spectrs	
#    import os
#    os.chdir(r'./exp_data')
#    
#    side_coefs = get_calibration_coefficients('clbr_coef_side.txt')
#    frame = read_files('tns.35-tsn.61',energy_scale = False, strip_convert = True)
#    ch_spectrs, en_spectrs = get_side_calibration_spectrs(frame, side_coefs)
#    visualize_side_spectr(1,ch_spectrs)
#    xmin,xmax = [],[]
#    for i in xrange(6):
#        xmin.append( int((5500 - side_coefs[1][i])/side_coefs[0][i]) )
#        xmax.append( int((10500 - side_coefs[1][i])/side_coefs[0][i]) )
#    xmin,xmax = np.array(xmin), np.array(xmax)
#                                    
#    Alpha_energies = [6040,6143.8,6264,6899.2,7137,7922,8699,9261]#[7922,8699,9261]
#    ind = 2
#    xpeaks,coef=calibration.calibrate_area(ch_spectrs[ind][1][1:],ch_spectrs[ind][0],xmin[ind],xmax[ind],threshold=0.08,sigma=3,visualize=True,energies=Alpha_energies)
#    print calibration.make_report(xpeaks,coef,i,energies=Alpha_energies)	

#   
##CALIBRATE FOCAL-SIDE ENERGY SPECTRS     
    class record: pass
    
    #energy scale check-calibration
    calibration_properties = record()
    calibration_properties.visualize=True
    calibration_properties.weighted_average_sigma = 14 #None#10
    calibration_properties.dlt = 170 # параметр для определения области соответствия расчетного положения пиков (соотв. пропорции калибровочных энергий) и реальных найденных пиков (метод spectrum_tools.search_peaks), если реальные пики не находятся вблизи расчетных 
                                    # в области dlt, то они заменяются на расчетные, в противном случае выбирается пик наиболее близкий к расчетному. see calibration.calibrate_area
    calibration_properties.energies = [6040,6143,6264,6899.2,7137,7922,8699,9261]#[7137,7922,8699,9261]#[6040,6143,6264,8699,9261]# 
    
    #calibration of aggregative focal-side spectr
    calibration_properties1 = record()
    calibration_properties1.visualize=True
    calibration_properties1.weighted_average_sigma = 10 #None#10
    calibration_properties1.dlt = 60 # параметр для определения области соответствия расчетного положения пиков (соотв. пропорции калибровочных энергий) и реальных найденных пиков (метод spectrum_tools.search_peaks), если реальные пики не находятся вблизи расчетных 
                                    # в области dlt, то они заменяются на расчетные, в противном случае выбирается пик наиболее близкий к расчетному. see calibration.calibrate_area
    calibration_properties1.energies = [6040,6143,6264,6899.2,7137,7922,8699,9261]#[7137,7922,8699,9261]#[6040,6143,6264,8699,9261]# 
    
    search_properties = record() #for smooth spectrum
    search_properties.widths=np.arange(1,10)
    search_properties.wavelet= 'ricker'
    search_properties.min_length=1.
    search_properties.min_snr=0.7
    search_properties.noise_perc= 0.3
    
    filter_properties = record()
    filter_properties.window_smooth=7
    filter_properties.smooth_wavelet='blackman'
    filter_properties.background_options='BACK1_ORDER8,BACK1_INCLUDE_COMPTON' 
    filter_properties.threshold = 4
    
    #калибровать отдельно по индексу
    def processing(ind_start,ind_stop,filename):
        for ind in xrange(ind_start,ind_stop):
            side_coefs[0,ind] = 2.4
            side_coefs[1,ind] = 0
            print 'strip: ',ind+1
            coefs,peaks = make_focal_side_calibration(frame,front_coefs,side_coefs,calibration_properties,calibration_properties1,search_properties,filter_properties,ind = ind)  
            print 'result: ',coefs
            try:
                a = coefs[1]
                b = coefs[0]
                print make_report(peaks,(a,b),ind+1,energies=calibration_properties.energies,filename = filename)[0]
            except:
                pass
            
    processing(1,2,filename = '/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/side_coef_clbr.txt')
        
        
#   CALIBRATE FOCAL FISSION SCALE
#    from read_files import read_files, get_front_spectrs
#    from calibrator import calibrator
#    import os
#    
#    #read data
#    os.chdir(r'/home/eastwood/codes/Python_Idle/data_processing/Files88_91')
#    frame = read_files('tsn.88-tsn.91',strip_convert=True)
#    hist,hsum  = get_front_spectrs(frame,tof=False,type_scale='Fchannel',visualize=True)
#    hist1,hsum1 = get_front_spectrs(frame,tof=True,type_scale='Fchannel',visualize=True)
#    
#    #select a spectrum including alpha-peaks from events without time-of-flight mark (tof=False, first 100 channels)
#    #and a peak of scattered ions which lies in area > 100 channel
#    ind = 7
#    def get_hist(i):
#        hist_probe = hist[i]
#        hist_probe[100:] = 0
#        hist_probe += hist1[i]
#        return hist_probe
#        
#    hist_probe = get_hist(ind)
#    
#    #create calibrator object to make calibrations by hand
#    Clbr = calibrator(energies = [6040,6143,6264,6899.2,7137,7922,8699,9261], calibration_function = exponential)
#        # while initiation of "calibrator" class one can set a list of energies and function for calibration of the final peak-energies list
#    Clbr.read(hist_probe.index,hist_probe)
#    Clbr.dl(25,500) #select an area from 25 to 500 channel
#    Clbr.set_energies([7922,8699,9261,196591])
#    
#    print(" Select 3 most right alpha-peaks (in area <100 channel) and high peak of scattered ions (nearby 450 channel) \n \
#and then put next commands into console panel to calculate calibration coefficients: \n \
#    Clbr.set_energies([7922,8699,9261,196591]) #energies in MeV of the peaks \n \
#    Clbr.calibrate(p0=(0,0.2,-10),filename='../Sep2014_calibrations/fission_clbr.txt',ind=ind)) #one could set initial coefficients for iterational calculating process and a filename for ouput report \n \
#if you have selected some needless excess peaks, you can delete them from the list using next commands: \n \
#    Clbr.show_points() # the command show a list of sorted selected peaks and their indexes \n \
#    Clbr.delete_points([5,6,7]) #set the indexes of excess peaks to delete them \n \
## the same works for energies, just use commands show_energies, delete_energies in the same way.\
#        ")
        
    
 
 
 