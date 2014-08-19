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



def calibrate_area(sample,xmin,xmax,search_peak_threshold=0.25,visualize=True):  
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
    xpeaks,ypeaks = search_peaks(xsample,hist,threshold=search_peak_threshold)
    #collecting the the largest 11 peaks
    dic_peaks = dict( zip(ypeaks,xpeaks) )
    dict1 = {}
    if len(dic_peaks)>=14: 
        count_peaks = 14
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
        
    #print xpeaks, ypeaks
    """selecting valid peaks"""   
    def get_numb(x,list1):
        for i,j in enumerate(list1):
            if j == x:
                return i
    
    if len(xpeaks) < 8:
        print 'Not enough peaks'
        print 'Peaks founded: ',xpeaks,ypeaks
        plt.show()
        raise ValueError('Not enough peaks')
        
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
    
    #print 'Peaks founded: ',xpeaks,ypeaks   
    if (xpeaks[-6]-xpeaks[-7])> 0.75*(xpeaks[-4]-xpeaks[-5]):
        #print 'YES'
        #print xpeaks
        xstart,xstop = get_numb(xpeaks[-7],xsample),get_numb(xpeaks[-6],xsample)
        #print xsample[xstart],xsample[xstop]
        #print xstart, xstop
        if xstop - xstart != 0:
            extra_xpeaks,extra_ypeaks = search_peaks(xsample[xstart:xstop],\
            hist[xstart:xstop],threshold = 0.3)
            #print extra_xpeaks
            #print '!!',extra_xpeaks,extra_ypeaks
            if len(extra_xpeaks)>=3:
                middle_peak_num = get_numb(sorted(extra_ypeaks)[-3], extra_ypeaks)
                xpeaks.insert(-6,extra_xpeaks[middle_peak_num])
                ypeaks.insert(-6,extra_ypeaks[middle_peak_num])
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
        
    def residuals(coef,y,x):
        return y - coef[0]*np.ones(len(x)) - coef[1]*x
   
    p0 = (0,2) #init coefficients    
    energies = np.array([6040,6143,6264,6899,7137,7922,8699,9265])
    xpeaks = np.array(xpeaks)
    solution = leastsq(residuals,p0, args = (energies,xpeaks) )[0]
    
    """OUTPUT"""
    #print 'Result (xpeaks,ypeaks):', xpeaks, ypeaks
    #print 'Solution:', solution
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
    

def calibrate_spectrum(filename,xmin,xmax,strips,output_file=None,**argv):
    """ 
    Function is for calibration of spectrs gotten from group of files
        filename - names of the files to open 'str'
        xmin,xmax - area of the valid spectrum int
        strips - numbers of strips to calibrate [list]
        argv - parameters to put-in inside  get_front_spectrs function {dict}
        output_file - name of ouput file contains a report 'str'
    Output: -
    """
    sample,sum_spectr = get_front_spectrs(filename,**argv)
    
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
            xpeaks,solution = calibrate_area(hist,xmin,xmax,visualize=True)
            print make_report(xpeaks,solution,ind,filename='clbr_coef.txt') 
        except ValueError:
            print '%d   Error occured: the spectr %d hasn\'t been calibrated \n'%(ind+1,ind+1)

def make_report(xpeaks,solution,ind,filename=None):
    energies = np.array([6040,6143,6264,6899,7137,7922,8699,9265]) 
    report = '\n%d    A = %2.3f ; B = %2.3f\n'%(ind+1,solution[0],solution[1])
    report += '%5s  %5s      %4s  %5s \n'%('Eexp','Ecal','Differ','Channel')
    for i,en in enumerate(energies):
        Ecalc = solution[0]+solution[1]*xpeaks[i]
        report += '%4.1f  %4.1f    %3.1f    %4.1f \n'%(en,Ecalc,Ecalc-en,xpeaks[i]) 
    
    if filename:
        f = open(filename,'a')
        f.write(report)
        f.close()
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
    """ INIT
    ind = 34
    sample,sum_spectr = get_front_spectrs('tsn.456-tsn.458',strip_convert=True,energy_scale=False,threshold=1,visualize = False)
    
    #PROCESSING
    xmin, xmax = 2000, 3700
    for ind in arange(0,47):#np.arange(len(sample)):#arange(0,3):
        hist = sample[ind]   
        

        #fig,ax = plt.subplots()
        #ax.plot(xsample,hist,linestyle = 'steps')
        #ax.plot(Bxp,Byp,'r',linestyle='steps')
        #ax.plot(Dxp,Dyp,'yD')
        #ax.plot([x_min,x_max],[30,30],'kD')
        #plt.show()
        # Собственно калибровка
        try:
            xpeaks,solution = calibrate_spectr(hist,xmin,xmax,visualize=True)
            print str(ind+1)+' '+make_report(xpeaks,solution) 
        except ValueError:
            print '%d   Error occured: the spectr %d hasn\'t been calibrated \n'%(ind+1,ind+1)
    """
        
    filename = 'tsn.456-tsn.458'
    arguments = {'strip_convert':'True','energy_scale':'False','threshold':1,'visualize':'False'}
    xmin, xmax = 2000, 3700
    strips = np.arange(0,47)
    calibrate_spectrum(filename,xmin,xmax,strips,**arguments)
    