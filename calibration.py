# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 15:08:07 2014

Code for processing automatic line-calibration of alpha-scale in our experiments.
Based on spectrum_tools.py 
 
@author: eastwoodknight
"""

from spectrum_tools import *
from read_files import * #get_front_spectrs
from scipy.optimize import leastsq
import numpy as np

Alpha_energies = [6040,6143.8,6264,6899.2,7137,7922,8699,9261]
 
def show_spectrs(xsample0,sample0,xsample,xmin,xmax,xpeaks,ypeaks,sample,solution,energies):
    fig,ax = plt.subplots(3,1,sharex=True)
    ax[0].set_title('Raw spectrum')
    ax[0].plot(xsample0,sample0,linestyle='steps') 
    ax[1].set_title('Processed spectrum with marked calibration peaks')
    ax[1].set_xlim([xmin,xmax])
    ax[0].set_ylim([0,max(ypeaks)*2.15])						
    ax[1].set_ylim([0,max(ypeaks)*1.25])
    ax[1].plot(xsample,sample,linestyle='steps',color='b',linewidth=3) 
    ax[1].plot(xpeaks,ypeaks,'ro',linewidth=4)
    ax[2].set_xlim([xmin,xmax])
    ax[2].set_title('Calibration function')
    ax[2].plot(xpeaks,energies,'ro',linewidth=4,label='Data points')
    ax[2].plot(xpeaks,solution[0]*np.ones(len(xpeaks))+solution[1]*xpeaks,linewidth=2,label='Fitting line')
    #print solution[0]*np.ones(len(xpeaks))+solution[1]*xpeaks
    ax[2].legend(loc='lower right')
    plt.show()
    
def filter_spectr(xmin,xmax,xsample,sample,properties): #window_smooth=7,smooth_wavelet='hanning',background_options='BACK1_ORDER4'):
    sample = np.array(sample)
    x_ind = (xsample >= xmin) & (xsample <= xmax)
    sample = sample[x_ind]
    xsample = (xsample[x_ind]).tolist() 
    if (sample == 0).sum() == len(sample):
        raise  Exception('No data in sample[xmin:xmax]')
    if properties.window_smooth: 
        sample = smooth(sample,properties.window_smooth,properties.smooth_wavelet)
        sample_background = background(sample,parameters=properties.background_options) #,BACK1_INCLUDE_COMPTON
        sample -= sample_background
        sample[ sample<0 ] = 0  
    return xsample,sample


 
        
def calibrate_area(xsample,sample,xmin,xmax,calibration_properties,filter_properties): #threshold=0.25,sigma=3,visualize=True,energies=Alpha_energies):  
    """ 
    Example of input parameter's classes:
    
    class record: pass
    calibration_properties = record()
    calibration_properties.threshold=0.12
    calibration_properties.sigma=2
    calibration_properties.visualize=True 
    calibration_properties.weighted_average_sigma = None
    calibration_properties.dlt = 30
    calibration_properties.energies = [6040,6143,6264,6899.2,7137,7922,8699,9261]
    filter_properties = record()
    filter_properties.window_smooth=3
    filter_properties.smooth_wavelet='hanning'
    filter_properties.background_options='BACK1_ORDER4'
    
    Return xpeaks, solution (result of linear fitting of spectrum)
    """
    xsample0 = xsample #make copies to visualize it later
    sample0 = sample
    
    #initiatian of properties
    try:
        sigma = calibration_properties.sigma
        threshold = calibration_properties.threshold
        visualize = calibration_properties.visualize
        energies = calibration_properties.energies
    except:
        raise Exception('Wrong calibration_properties object')
    
    #filter data - smooth and delete a background
    if filter_properties:
        xsample,sample = filter_spectr(xmin,xmax,xsample,sample,filter_properties)#5,smooth_wavelet='hanning',background_options='BACK1_ORDER2')
    
    #find peaks
    xpeaks,ypeaks = search_peaks(xsample,sample,sigma=sigma,threshold=threshold) 
    
    #delete too close peaks and sort them
    indx = np.concatenate((np.abs(np.diff(xpeaks))>5 ,[True]))  
    xpeaks,ypeaks= xpeaks[indx],ypeaks[indx]
    ypeaks = ypeaks[xpeaks.argsort()] #sort lists of peaks
    xpeaks = sorted(xpeaks)
    xpeaks,ypeaks = np.array(xpeaks),np.array(ypeaks)
    print xpeaks
    
    #selecting valid peaks
    spectr_peak_dists = np.diff(np.array(energies,dtype = np.float64) ) / (energies[-1]-energies[0])
    spectr_length = (xpeaks[-1] - xpeaks[-2])/spectr_peak_dists[-1] 
    spectr_peak_dists1 = spectr_peak_dists*spectr_length # distances from right edge to points, proportion between them is the same as for energy calibration spectr
    spectr_peak_dists1 = xpeaks[-1] - spectr_peak_dists1[::-1].cumsum()
    x,y = [],[]
    x.append(xpeaks[-1])
    y.append(ypeaks[-1])
    dlt = calibration_properties.dlt    
    def find_closest(msv,m0,dlt): #find a point from msv array, which is the closest to m0 inside dlt diapason
        msv=abs(msv-m0)
        i = msv.argmin()
        if msv[i] < dlt:
            return i
        else:
            return False
                
    for i in xrange(len(spectr_peak_dists1)):
        l =spectr_peak_dists1[i]
        peak_ind = find_closest(xpeaks,l,dlt)
        if not peak_ind:
            ind = (xsample > l).argmax()
            x.append(xsample[ind])
            y.append(sample[ind])
            continue
        x.append(xpeaks[peak_ind])  
        y.append(ypeaks[peak_ind]) 
        spectr_length += l - x[-1] 
        spectr_peak_dists1 = spectr_peak_dists*spectr_length
        spectr_peak_dists1 = xpeaks[-1] - spectr_peak_dists1[::-1].cumsum()
               
    x.reverse()
    y.reverse()
    xpeaks,ypeaks= np.array(x,dtype=np.float64),np.array(y,dtype=np.float64)  
    if len(xpeaks) < len(energies):
        print 'Not enough peaks'
        print 'Peaks founded: ',xpeaks,ypeaks
        plt.show()
        raise ValueError('Not enough peaks') 
         
    #correct xpeaks by calculating a weighted average of x using +-2*sigma window   
    if calibration_properties.weighted_average_sigma:
        k = calibration_properties.weighted_average_sigma
        print 'start wa',k
        xpeaks1 = []
        for x in xpeaks:
            #print xsample, x in xsample			
            ind = xsample.index(x)  
            hist_sum = sample[ind-k:ind+k+1].sum()
            xpeaks1.append( (sample[ind-k:ind+k+1]*xsample[ind-k:ind+k+1]).sum()/hist_sum )
        xpeaks = np.array(xpeaks1) 
        
    #fitting by line
    def residuals(coef,y,x):
        return y - coef[0]*np.ones(len(x)) - coef[1]*x
    p0 = (0,2) #init coefficients    
    
    xpeaks,ypeaks= np.array(xpeaks,dtype=np.float64),np.array(ypeaks,dtype=np.float64)
    solution = leastsq(residuals,p0, args = (energies,xpeaks) )[0]
    
    #output
    if visualize:
        show_spectrs(xsample0,sample0,xsample,xmin,xmax,xpeaks,ypeaks,sample,solution,energies)
    print 'pr',xpeaks
    return xpeaks,solution    
    
    
    
def calibrate_spectrum(hist,xmin,xmax,strip_ind='[1,2,3]',output_file=None,threshold=0.25,sigma=3,visualize=True,energies=Alpha_energies):		
	exec('strips = hist.columns['+strip_ind+']')
	for i in hist[strips]:
	    xeaks,coef=calibrate_area1(hist.index,np.array(hist[i]),xmin,xmax,threshold=threshold,sigma=sigma,visualize=visualize,energies=energies)
	    print make_report(xpeaks,coef,i,filename=output_file,energies=Alpha_energies)				
				
"""    
def calibrate_spectrum(filename,xmin,xmax,strips,output_file=None,args={},search_agrs={}):
     
    Function is for calibration of spectrs gotten from group of files
        filename - names of the files to open 'str'
        xmin,xmax - area of the valid spectrum int
        strips - numbers of strips to calibrate [list]
        argv - parameters to put-in inside  get_front_spectrs function {dict}
        output_file - name of ouput file contains a report 'str'
    Output: -
    
    sample,sum_spectr = get_front_spectrs(filename,**args)
    
    if output_file:
        filename = output_file
    else:
        filename ='clbr_coef_Sep2014.txt'
        
    for ind in strips:#np.arange(len(sample)):#arange(0,3):
        hist = sample[ind]   
   # Визуализация     
        #fig,ax = plt.subplots()
        #ax.plot(xsample,hist,linestyle = 'steps')
        #ax.plot(Bxp,Byp,'r',linestyle='steps')
        #ax.plot(Dxp,Dyp,'yD')
        #ax.plot([x_min,x_max],[30,30],'kD')
        #plt.show()
    # Собственно калибровка
        try:
            xpeaks,solution = calibrate_area(hist,xmin,xmax,**search_agrs)
            print make_report(xpeaks,solution,ind,filename) 
        except ValueError:
            print '%d   Error occured: the spectr %d hasn\'t been calibrated \n'%(ind+1,ind+1)
        except IndexError:
            print '%d   Error occured: possibly wrong area of calibration (xmin,xmax) \n'%(ind+1)
        except KeyError:
            print 'No index: %d' %(ind+1)
"""

def make_report(xpeaks,solution,ind,filename=None,energies=Alpha_energies):
    #energies = np.array([6040,6143,6264,6899,7137,7922,8699,9265]) 
    report = '\n%d    A = %2.5f ; B = %2.1f\n'%(ind,solution[1],solution[0])
    report += '%5s  %5s      %4s   %5s \n'%('Eexp','Ecal','Differ','Channel')
    S = 0
    for i,en in enumerate(energies):
        Ecalc = solution[0]+solution[1]*xpeaks[i]
        report += '%4.1f  %4.1f    %-5.1f    %-4.1f \n'%(en,Ecalc,Ecalc-en,xpeaks[i]) 
        S += (Ecalc-en)**2
    
    if filename:
        f = open(filename,'a')
        f.write(report)
        f.close()
    print 'S/n = ',S**0.5/len(energies)
    return report   




def find_calibration_area(hist):

    hist[hist > 300] = 300
    xsample = np.arange(len(hist))
    # Отбор высоких (Byp), "разграничивающих" пиков и мелких, которые позволяеют выявить (Dxp)наибольшую плотность скопления
    #подходящих пиков. Исп-ся пороги в соотношении примерно 3:1
    Bxp, Byp = search_peaks(xsample,hist,threshold=0.35)
    Dxp, Dyp = search_peaks(xsample,hist,threshold=0.07)
    # Выбор наиболее плотной области между разграничивающими пиками и определение границ этой области 
    #  в этой области и будет проводится автокалибровка
   
    Bxp = Bxp.tolist()
    Bxp.append(max(xsample))
    i = 0
    xstart_list =[]
    xstop_list =[]
    density_list =[]
    max_i = len(Dxp)-1

    x_border = 3700
    for xp in Bxp:
        density = 0
        j = i

        while Dxp[i] < xp: 
            i += 1
            density+=1
            if (i == max_i) or (Dxp[i] > x_border):
                break

        if i > j+2:   
            xstart_list.append(Dxp[j+1])
            xstop_list.append(Dxp[i])
            density_list.append(density)

        if i == max_i or Dxp[i] > x_border:
            break   

    #print 'Borders, density peaks: ', Bxp, Dxp
    #print 'xstart, xstop', xstart_list, xstop_list
    Bxp.pop()
    Bxp = np.array(Bxp)
    #print xstart_list, xstop_list, density_list
    numb = density_list.index(max(density_list))
    xmin,xmax = xstart_list[numb]+50, xstop_list[numb]+50
    
    return xmin,xmax   
        
        
        
if __name__ == '__main__':
#    
#    #INIT
#    ind = 34
#    sample,sum_spectr = get_front_spectrs('tsn.456-tsn.458',strip_convert=True,energy_scale=False,threshold=1,visualize = False)
#    
#    #PROCESSING
#    xmin, xmax = 2000, 3700
#    for ind in arange(1,47):#np.arange(len(sample)):#arange(0,3):
#        hist = sample[ind]   
#        
#
#        #fig,ax = plt.subplots()
#        #ax.plot(xsample,hist,linestyle = 'steps')
#        #ax.plot(Bxp,Byp,'r',linestyle='steps')
#        #ax.plot(Dxp,Dyp,'yD')
#        #ax.plot([x_min,x_max],[30,30],'kD')
#        #plt.show()
#        # Собственно калибровка
#        try:
#            xpeaks,solution = calibrate_spectr(hist,xmin,xmax,visualize=True)
#            print str(ind+1)+' '+make_report(xpeaks,solution) 
#        except ValueError:
#            print '%d   Error occured: the spectr %d hasn\'t been calibrated \n'%(ind+1,ind+1)
    
      
    filename = 'tsn.456-tsn.461'
    arguments = {'strip_convert':'True','energy_scale':'False','threshold':1,'visualize':'False'}
    search_arg = {'threshold':0.28,'visualize':'True'}
    xmin, xmax = 2050, 3700
    #strips = np.arange(1,4)
    #calibrate_spectrum(filename,xmin,xmax,strips,output_file = 'clbr_Shumeiko_2-Oct.txt',args=arguments,search_agrs=search_arg)#,output_file = 'clbr_Shumeiko_2-Oct.txt')
    
    hist,sum1 = get_front_spectrs(filename,strip_convert=True,threshold=0.04)
    x,y = hist.shape
    x = np.arange(1,x+1)
    xpeaks,coef=calibrate_area(x,np.array(hist[2]),xmin,xmax,threshold=0.15,sigma=5,visualize=True,energies=Alpha_energies)
    print make_report(xpeaks,coef,2,energies=Alpha_energies)
    
    Alpha_energies = [6040,6143.8,6264,6899.2,7137,7922,8699,9261] #'clbr_coef_back_21Oct.txt',
    #Alpha_energies = [7137,7922,8699,9261]
    #hist,sum1 = get_front_spectrs(filename,strip_convert=True,window=False,threshold=0.04,visualize=False)    
    #calibrate_spectrum(hist,680,1220,strip_ind='0:3',output_file=None,threshold=0.15,sigma=3,visualize=True,energies=Alpha_energies)
	