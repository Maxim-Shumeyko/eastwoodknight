# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 20:45:04 2015

@author: eastwood
"""

from read_files import *
from read_american_format_cy import *
from calibration import *
 
def visualize_side_spectr(strip,xsample,spectrs):
    sample1 = smooth(spectrs[strip],7,'hanning')
    sample1[sample1<0] = 0
    sample = sample1 - background(sample1,niter=3,parameters='BACK1_ORDER8,BACK1_INCLUDE_COMPTON')
    sample[sample<0] = 0
    fig,ax = plt.subplots(3,1,sharex=True,sharey=True)		
    ax[0].plot(xsample,sample,linestyle='steps')
    ax[1].plot(xsample,sample1,linestyle='steps')
    ax[2].plot(xsample,spectrs[strip],linestyle='steps')
    plt.show()
    return sample
 
def calibrate_side_calibration_spectrs(frame,side_coefs,front_coefs,calibration_properties,filter_properties,step=10):
    xsample,ch_spectrs, en_spectrs = get_side_calibration_spectrs(frame, side_coefs,front_coefs,step=step)
    
    #calculate calibration boundaries of spectrums
    xmin,xmax = [],[]
    for i in xrange(6):
        xmin.append( int((5500 - side_coefs[1][i])/side_coefs[0][i]) )
        xmax.append( int((10500 - side_coefs[1][i])/side_coefs[0][i]) )
    xmin,xmax = np.array(xmin), np.array(xmax)
            
    ind = 2 #select a strip to visualize 
    if calibration_properties.visualize:
        #ch_spectrs[ind] = 
        visualize_side_spectr(ind,xsample,en_spectrs)
    xpeaks,coef=calibrate_area(xsample,ch_spectrs[ind],xmin[ind],xmax[ind],calibration_properties,filter_properties)
    print make_report(xpeaks,coef,ind,energies=calibration_properties.energies)
    return coef
 
 
if __name__ == '__main__':
    
    
#CONSTRUCT AND SHOW FOCAL-SIDE ENERGY SPECTRS 
    
    #go to the data folder
    import os
    os.chdir(r'./exp_data')
        
    side_coefs = get_calibration_coefficients('clbr_coef_side.txt')
    front_coefs = get_calibration_coefficients('clbr_coef_front.txt')
    frame = read_files('tns.35-tsn.50',energy_scale = False, strip_convert = True) #-tsn.61
    side_coefs[0,2] = 2.5
    side_coefs[1,2] = 0
    
    #construct calibration objects of properties
    class record: pass
    calibration_properties = record()
    calibration_properties.threshold=0.2
    calibration_properties.sigma=2
    calibration_properties.visualize=False
    calibration_properties.weighted_average_sigma = 4
    calibration_properties.dlt = 30 # параметр для определения области соответствия расчетного положения пиков (соотв. пропорции калибровочных энергий) и реальных найденных пиков (метод spectrum_tools.search_peaks), если реальные пики не находятся вблизи расчетных 
                                    # в области dlt, то они заменяются на расчетные, в противном случае выбирается пик наиболее близкий к расчетному. see calibration.calibrate_area
    calibration_properties.energies = [7137,7922,8699,9261]#[6040,6143,6264,6899.2,7137,7922,8699,9261]#[6040,6143,6264,8699,9261]# 

    filter_properties = record()
    filter_properties.window_smooth=7
    filter_properties.smooth_wavelet='hanning'
    filter_properties.background_options='BACK1_ORDER2' 
    
    for i in xrange(9):
        coef = calibrate_side_calibration_spectrs(frame,side_coefs,front_coefs,calibration_properties,filter_properties,step=10)
        side_coefs[0,2] = coef[1]
        side_coefs[1,2] = coef[0]
    
 
 