# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 20:45:04 2015

@author: eastwood
"""

from read_files import *
from read_american_format_cy import *
from calibration import *
import numpy as np
from time import time
 
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
    a = side_coefs[0,ind-1]*coef[1] #A coef[1]; clbr_func = (A*x + B)
    b = side_coefs[1,ind-1]*coef[1] + coef[0] #B coef[0] 
    side_coefs[0,ind-1] = a
    side_coefs[1,ind-1] = b
    control_S2 = sum((xpeaks - np.array(calibration_properties.energies))**2)**0.5/len(xpeaks)
    print ' Control: ', control_S2
    print
    if control_S2 > 350:
        raise Bad_result, 'wrong coefficients or bad calibration_properties'
    #calculate calibration boundaries of spectrums
    xmin = int((4500 - b)/a) # a == side_coefs[0][i] 
    xmax = int((9600 - b)/a)  # b == side_coefs[1][i] 
    xsample,ch_spectrs, en_spectrs = get_side_calibration_spectrs(frame, side_coefs,front_coefs,step=step,strip=ind)
    xpeaks,coef=calibrate_area(xsample,ch_spectrs,xmin,xmax,calibration_properties1,filter_properties,search_properties) 
    return coef,xpeaks, control_S2
 
   
#CONSTRUCT AND SHOW FOCAL-SIDE ENERGY SPECTRS 
    
def make_focal_side_calibration(frame,front_coefs,side_coefs,calibration_properties,calibration_properties1,search_properties,filter_properties,ind = 1,iterations=3):   
    """
    Function for making autoumatic calibration of focal-side spectrums strip by strip.
    Output: coefs (A,B), list of peaks (which were used in calibration)
    """
    
    control, coefs, peaks = [],[],[]    
    side_coefs_start = np.array(side_coefs)
    coef = [side_coefs[0,ind-1],side_coefs[1,ind-1] ]
    xpeaks = False
    for i in xrange(iterations):
        try:
            coefs.append(coef)
            peaks.append(xpeaks)
            coef,xpeaks,control_S2 = calibrate_side_calibration_spectrs(frame,side_coefs,front_coefs,calibration_properties,calibration_properties1,filter_properties,search_properties,ind,flag=True,step=10)
            control.append(control_S2)
            side_coefs[0,ind-1] = coef[1] #A coef[1]
            side_coefs[1,ind-1] = coef[0] #B coef[0] 
        except (Calibration_Error,Bad_result),e:
            while len(coefs) > len(control):
                coefs.pop()
            side_coefs[0,ind-1] += 0.05*np.random.randn()
            side_coefs[1,ind-1] += 20*np.random.randn()
        print 'coefs:',side_coefs[0,ind-1],side_coefs[1,ind-1]
    if len(coefs) == 0:
        print 'Bad initial coefficients'
        return
    ind1 = np.argmin(np.array(control))   
    side_coefs[0,ind-1] = coefs[ind1][1] #A coef[1]
    side_coefs[1,ind-1] = coefs[ind1][0] #B coef[0]
    print 'Control: ', control[ind1]
    return (coefs[ind1][1], coefs[ind1][0]), peaks[ind1]
    
if __name__ == '__main__':

# READ DATA
#go to the data folder and read kamak data
#    import os

#    side_coefs = get_calibration_coefficients('clbr_coef_side.txt')
#    print 'right coefs:',side_coefs,'\n'
#    front_coefs = get_calibration_coefficients('clbr_coef_front.txt')
#    os.chdir(r'./exp_data1')
#    frame = read_files('tsn.35-tsn.61',energy_scale = True, strip_convert = True)#\
#        clbr_front_filename='/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/alpha_clbr.txt') #-tsn.61
#READ DATA AND PLOT DISTRIBUTIONS 
#    os.chdir('May_2015')
#    frame1 = read_files('tsn.6-tsn.7',strip_convert=True)
#    hist,hsum = get_front_spectrs(frame1,threshold=0.04)
#    hist_b,hsum_b = get_back_spectrs(frame1,threshold=0.04)
#    hist_side,hsum_side = get_front_spectrs(frame1,id_mark='==3',threshold=0.04)
#READ AMERICAN DATA AND PLOT DISTRIBUTIONS 
#     # example of counts and back distribution for a large file
#    os.chdir('May_2015_AM')
#    frames = read_am.OpenAmFile('24May_15_04_cal.bin',chunk_size=1000000).get_frames()
#    frame = frames.next()
#    hist,hist_s = get_am_back_spectrs(frame,threshold = 0.03,visualize=False)
#    counts= frame.groupby('back_strip_even').count()['strip']
#    for frame in frames:
#        counts = frame.groupby('back_strip_even').count()['strip']
#        
#        hist1,hist_s1 = get_am_back_spectrs(frame,threshold = 0.03,visualize=False)
#        hist += hist1
#        hist_s += hist_s1
#    fig,ax = plt.subplots()
#    ax.plot(counts.index,counts,linestyle='steps')
#    visualize_spectrum(hist,hist_s)
#    Simple example
#    frame1 = read_am.OpenAmFile('22May_cal.bin',chunk_size=10000000).get_frames().next()
#    hist,hsum = get_front_spectrs(frame1,xsize=16384,threshold=0.05)
#    hist_b,hsum_b = get_am_back_spectrs(frame1,threshold=0.05)
#    hist_side,hsum_side = get_front_spectrs(frame1,id_mark='==3',xsize=16384,threshold=0.05)
#MAKE DISTRIBUTIONS and CALIBRATE ALPHA-SCALE OF AMERICAN FILES
#    files_list = ['22May_15_03_cal1.bin','22May_15_04_cal.bin']    
#    os.chdir('May_2015_AM')
#    final_hist,final_hist_s = pd.DataFrame(np.zeros((16382, 128)), columns=np.arange(1,129)),np.zeros(16382)
#    for file_ in files_list:
#        
#        print file_,'is processing'
#        time_start = time()
#        frames = read_am.OpenAmFile(file_,chunk_size=1000000).get_frames(time_stop='2015-05-22 17:52')
#        frame = frames.next()
#        hist,hist_s = get_am_back_spectrs(frame,threshold = 1.,visualize=False)
#        for column in hist.columns:
#            final_hist[column] += hist[column]
#        final_hist_s += hist_s
#        break
#        #counts= frame.groupby('back_strip_even').count()['strip']
#        for frame in frames:
#            #counts = frame.groupby('back_strip_even').count()['strip']          
#            hist1,hist_s1 = get_am_back_spectrs(frame,threshold = 1.,visualize=False)
#            hist += hist1
#            hist_s += hist_s1 
#            del hist1,hist_s1
#            #del frame
#        print 'shape',hist.shape, final_hist.shape
#        for column in hist.columns:
#            final_hist[column] += hist[column]
#        final_hist_s += hist_s
#        del hist,hist_s
#        print 'time of processing:',time()-time_start
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
#FRONT SPECTRS
#    hist,sum_spectr = get_front_spectrs(frame,energy_scale=True) #'energy_scale' provides energy spectr 
#    fig,ax = plt.subplots()
#    ax.plot(hist.index,sum_spectr,linestyle='steps')
#    plt.show()
#CALIBRATE FRONT (OR BACK) ALPHA SPECTRS 
#    import os
#    os.chdir(r'./May_2015')
#    #read data 
#    filenames = 'tsn.6,tsn.7,tsn.8,tsn.9,tsn.10,tsn.11,tsn.12'.split(',')
#    start = time()
#    #frame = read_files(filenames,strip_convert=False)
#    print time()-start
#    start = time()
#    hist,hist_b = pd.DataFrame(np.zeros((48,8192)),index=np.arange(1,49)),pd.DataFrame(np.zeros((128,8192)),index=np.arange(1,129))
#    hsum,hsum_b = np.zeros(8192),np.zeros(8192)
#    for file_ in filenames:
#        hist1,hsum1 = get_front_spectrs(file_,threshold=0.04,visualize=False)
#        hist += hist1
#        hsum += hsum1
#        hist_b1,hsum_b1 = get_back_spectrs(file_,threshold=0.04,visualize=False)
#        hist_b += hist_b1
#        hsum_b += hsum_b1
#    filenames = 'tsn.6-tsn.12'
#    hist_b,hsum_b= get_back_spectrs(filenames,threshold=0.04)
#    print time()-start
#    xmin, xmax = 1100, 2100
#    def get_calibration_properties():
#        #set calibration parameters
#        class record: pass
#        calibration_properties = record()
#        calibration_properties.visualize= True
#        calibration_properties.weighted_average_sigma = 5
#        calibration_properties.dlt = 12 # параметр для определения области соответствия расчетного положения пиков (соотв. пропорции калибровочных энергий) и реальных найденных пиков (метод spectrum_tools.search_peaks), если реальные пики не находятся вблизи расчетных 
#                                        # в области dlt, то они заменяются на расчетные, в противном случае выбирается пик наиболее близкий к расчетному. see calibration.calibrate_area
#        calibration_properties.energies = [6040,6143,6264,6899.2,7137,7922,8699,9261]#[7137,7922,8699,9261]#[6040,6143,6264,8699,9261]# 
#    
#        search_properties = record() #for noisy spectrum
#        search_properties.widths=np.arange(1,5)
#        search_properties.wavelet= 'ricker'
#        search_properties.min_length=1.
#        search_properties.min_snr=0.6
#        search_properties.noise_perc=0.2 
#        
#        filter_properties = record()
#        filter_properties.window_smooth=7
#        filter_properties.smooth_wavelet='blackman'
#        filter_properties.background_options='BACK1_ORDER8,BACK1_INCLUDE_COMPTON' 
#        filter_properties.threshold= 8
#        return calibration_properties, search_properties, filter_properties
#        
#    calibration_properties, search_properties, filter_properties = get_calibration_properties()
##    #choose strips and calibrate them
#    start,stop = 1,10
#    #output_filename = '/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/alpha_back_clbr.txt'
#    output_filename = None
#    xpeaks, coefs = calibrate_spectr(start,stop,xmin,xmax,hist_b,calibration_properties,filter_properties,search_properties,output_filename = output_filename)
## write coefs to file
#    f = open('coef_alpha.coef','a')
#    for i in coef:
#         f.write('%2.4f %4.1f\n' %(i[1],i[0])) 
#    f.close()    
				
##GET AND SHOW FOCAL-SIDE ENERGY SPECTRS
#    import os
#    os.chdir(r'./exp_data1')
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
##CALIBRATE ALPHA FOCAL-SIDE ENERGY SPECTRS   

    os.chdir(r'./exp_data1')
    frame = read_files('tsn.35-tsn.61',energy_scale = False, strip_convert = True)
    side_coefs = get_calibration_coefficients('/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/side_coef_clbr.txt')
    front_coefs = get_calibration_coefficients('/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/alpha_clbr.txt')
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
    def processing(frame,ind_start,ind_stop,filename,side_coefs=side_coefs):
        for ind in xrange(ind_start,ind_stop):
            side_coefs[0,ind-1] = 2.4
            side_coefs[1,ind-1] = 0
            print side_coefs[0,ind-1],side_coefs[1,ind-1]
            print 'strip: ',ind
            coefs,peaks = make_focal_side_calibration(frame,front_coefs,side_coefs,calibration_properties,calibration_properties1,search_properties,filter_properties,ind = ind)  
            print 'result: ',coefs
            try:
                a = coefs[1]
                b = coefs[0]
                print make_report(peaks,(a,b),ind,energies=calibration_properties.energies,filename = filename)[0]
            except:
                pass
            
    processing(frame,1,2,filename = '/home/eastwood/codes/Python_Idle/data_processing/May_2015/side_coef_clbr.txt')
#        
#        
##CALIBRATE FOCAL FISSION SCALE
#    from read_files import read_files, get_front_spectrs
#    from calibrator import calibrator,exponential
#    import os
#    
#    #read data
#    os.chdir(r'/home/eastwood/codes/Python_Idle/data_processing/Files88_91')
#    frame = read_files('tsn.88-tsn.91',strip_convert=True)
#    hist,hsum  = get_front_spectrs(frame,tof=False,type_scale='fission_channel',visualize=True)
#    hist1,hsum1 = get_front_spectrs(frame,tof=True,type_scale='fission_channel',visualize=True)
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
        
        
#CALIBRATE BACK FISSION SCALE
#    from read_files import read_files, get_back_spectrs
#    from calibrator import calibrator,exponential
#    import os
#    
#    #read data
#    os.chdir(r'/home/eastwood/codes/Python_Idle/data_processing/Files88_91')
#    frame = read_files('tsn.88-tsn.91',strip_convert=True)
#    hist,hsum  = get_back_spectrs(frame,tof=False,fission_scale=True,visualize=True)
#    hist1,hsum1 = get_back_spectrs(frame,tof=True,fission_scale=True,threshold=0.4,visualize=True)
#    
#    #select a spectrum including alpha-peaks from events without time-of-flight mark (tof=False, first 100 channels)
#    #and a peak of scattered ions which lies in area > 100 channel
#    ind = 1
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
#    Clbr.dl(25,700) #select an area from 25 to 500 channel
#    Clbr.set_energies([7922,8699,9261,196591])
#    
#    print(" Select 3 most right alpha-peaks (in area <100 channel) and high peak of scattered ions (nearby 450 channel) \n \
#and then put next commands into console panel to calculate calibration coefficients: \n \
#    Clbr.set_energies([7922,8699,9261,196591]) #energies in MeV of the peaks \n \
#    Clbr.calibrate(p0=(0,0.2,-10),filename='../Sep2014_calibrations/fission_back_clbr.txt',ind=ind) #one could set initial coefficients for iterational calculating process and a filename for ouput report \n \
#if you have selected some needless excess peaks, you can delete them from the list using next commands: \n \
#    Clbr.show_points() # the command show a list of sorted selected peaks and their indexes \n \
#    Clbr.delete_points([5,6,7]) #set the indexes of excess peaks to delete them \n \
## the same works for energies, just use commands show_energies, delete_energies in the same way.\
#        ")
        
        
#  DISTRIBUTIONS ON FISSION SCALE
#    os.chdir(r'./Files88_91')
#    os.chdir(r'./exp_data1')
#    frame = read_files('tsn.35-tsn.61',energy_scale = True, \
#    #coefs,strips1,Fchannel1,res = read_file('tsn.35',energy_scale = True,\
#        fission_energy_scale = True,\
#        clbr_front_filename='/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/alpha_clbr.txt',\
#        clbr_fission_side_filename='/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/fission_side_coef_clbr1.txt',\
#        clbr_fission_back_filename='/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/fission_back_clbr.txt',\
#        clbr_fission_front_filename='/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/fission_clbr.txt')
#    hist,hsum  = get_fission_spectrs(frame,xsize = 20000, window = 50)
##    hist_side,hsum_side = get_fission_spectrs(frame,id_mark='side')
##    hist,hsum = get_back_spectrs(frame,fission_scale=True)
##  FOCAL-SIDE FISSION SPECTRS
##    os.chdir(r'./Files88_91')
##    os.chdir(r'./exp_data1')
#    frame = read_files('tsn.35-tsn.61',strip_convert = True)
#    side_coefs = get_calibration_coefficients('/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/fission_side_coef_clbr1.txt')
#    front_coefs = get_fission_calibration_coefficients('/home/eastwood/codes/Python_Idle/data_processing/Sep2014_calibrations/fission_clbr.txt')
##
##    side_coefs[0][0] = 12.
##    side_coefs[1][0] = 0.
#    
#    ind = 2
#    xsample,ch_spectrs, xsample_en,en_spectrs = get_side_calibration_spectrs(frame, side_coefs,front_coefs,fission_scale=True,strip=ind)
#    fig,ax = plt.subplots()
#    ax.plot(xsample_en,en_spectrs,linestyle='steps')
#    plt.show()
#    
#    ind = 3
#    xsample,ch_spectrs, xsample_en,en_spectrs = get_side_calibration_spectrs(frame, side_coefs,front_coefs,fission_scale=True,strip=ind)
#    fig,ax = plt.subplots()
#    ax.plot(xsample_en,en_spectrs,linestyle='steps')
#    plt.show()
      