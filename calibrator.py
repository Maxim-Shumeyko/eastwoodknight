# -*- coding: utf-8 -*-
"""
Created on Wed May 14 14:40:14 2014
This code is for simple calibration of gaussian peaks by hand.
@author: eastwoodknight
"""
import numpy as np
from mpl_cmd import fit_commander
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import leastsq

#fitting method function

def line(p,x):
    return p[0]*np.ones(len(x)) + p[1]*x 

def exponential(p,x):
    return p[0]*np.ones(len(x)) + p[1]*np.exp(-x/p[2]) 

def gaussian(p,x):
    return p[0]*np.exp( -((x - p[1])/p[2])**2 /2 )/np.sqrt( 1 / np.pi) #/ p[2]

        
def fit_method(x,y,func=gaussian):
    
    
    def residuals(p,y,x):
        return y - func(p,x)
    
#    #derivative
#    def Jac(p0,y,x):
#        dfdA = np.exp( -((x - p0[1])/p0[2])**2 /2 )/np.sqrt( 1 / np.pi) / p0[2]
#        dfdx0 = -( x - p0[1] )*func(p0,x)/2/p0[2]**2 #* residuals(p0,y,x) 
#        dfdsigm = (1/p0[2] - ( x - p0[0])**2 /(p0[2]**3) /2 )*func(p0,x)#* residuals(p0,y,x) 
#        deriv = np.c_[dfdA,dfdx0, dfdsigm]
#        return deriv
      
    #initial patameters x0, sigm, A(mplitude)
    p0 = []
    x0 = (x*y).sum()/y.sum()
    sigm0 = (( x - x0*np.ones(len(x)) )**2).sum()/ (len(x) - 2)
    A = y.sum() 
    p0.append(A)
    p0.append(x0)
    p0.append(sigm0)
    p0 = tuple(p0)
    
    good_calibration = True   
    try:
        solution = leastsq(residuals, p0,# Dfun = Jac,
                args = (y,x) )[0] #A,x0,sigm
        A,x0,sigm = solution  
    except Exception,e:
        print e
        A,x0,sigm = None,None,None
    finally:
        if (x0<min(x))or(x0>max(x)):
            print 'bad dataset for gaussian fitting'
            x0 = (x*y).sum()/y.sum()
            good_calibration = False
    return A,x0,sigm,good_calibration


def make_report(xpeaks,solution,ax,calibration_function = line,filename=None,energies=None,ind=1):
    #energies = np.array([6040,6143,6264,6899,7137,7922,8699,9265]) 
    report = '\n Solution:'
    
    if calibration_function == line:
        report = '\n%d    A = %2.5f ; B = %2.1f\n'%(ind,solution[1],solution[0])
    else:
        report = '\n%d    ' %(ind)
        msv1 = ['A','B','C','D','E','F']
        for i,coef in enumerate(solution):
            report += str(msv1[i])+' = '+str(coef)+' ; '
        
    report += '\n%5s  %5s      %4s   %5s \n'%('Eexp','Ecal','Differ','Channel')
    S = 0
    if calibration_function == line:
        func = lambda p,x: p[0] + p[1]*x
    elif calibration_function == exponential:
        func = lambda p,x: p[0] + p[1]*np.exp(-x/p[2]) 
    else:
        print 'Calibration function is unsuitable'
        return
     
    ypeaks = []       
    for i,en in enumerate(energies):
        Ecalc = func(solution,xpeaks[i])#solution[0]+solution[1]*xpeaks[i]
        report += '%4.1f  %4.1f    %-5.1f    %-4.1f \n'%(en,Ecalc,Ecalc-en,xpeaks[i]) 
        S += (Ecalc-en)**2
    
    xpeaks = np.linspace(min(xpeaks),max(xpeaks),500)
    ax.plot(xpeaks,func(solution,xpeaks),'g--')
    plt.show()
    if filename:
        f = open(filename,'a')
        f.write(report)
        f.close()
    report += 'S/n = %3.1f \n' % (S**0.5/len(energies))
    return report#, (S**0.5/len(energies))   
    
#main class     
class calibrator:
    
    def __init__(self,x = None, y = None,energies=None,calibration_function = line):
        
        fig, ax = plt.subplots()
        self.__fig, self.__ax = fig,ax
        self.read(x,y)
        self.__calibration_x = []
        self.__calibration_y = []
        self.calibration_function = calibration_function
        if energies:
            self.__energies = energies
        else:
            self.__energies = []
        #processing method
        def method(x,y):
            A,x0,sigm,good_calibration = fit_method(x,y)
            func = gaussian
            #delete previous calibrated peak from figure
            if len(self.__ax.lines) > 1:
                line = self.__ax.lines.pop(-1)
                del line
            #check if calibration was good    
            if good_calibration:
                print (A,x0,sigm)
                y0 = func((A,x0,sigm),x0)
                self.__ax.plot(x, func((A,x0,sigm),x), 'ro-')
            else:
                print x0
                y0 = y[np.argmax(x>=x0)]
                self.__ax.plot(x0, y0, 'ro-')
                
            self.__calibration_x.append(x0)
            self.__calibration_y.append(y0)   
            
            
        #connect event processor to the axes (ploting window)
        br = fit_commander(self.__ax,method)
        self.__fig.canvas.mpl_connect('button_press_event', lambda event:br.onclick(event))
    
    
    def set_energies(self,energies):
        if not energies:
            print 'Empty list'
        else:
            self.__energies = energies
      
      
    def delete_points(self,ind_list): 
        """
        delete points from calibration list
        set indexes of peaks to delete or empty list [] to delete them all
        """
        if ind_list == []:
            ind_list = range(len(self.__calibration_x))
        self.__calibration_x = filter(lambda x: self.__calibration_x.index(x) not in ind_list,self.__calibration_x)
        self.__calibration_y = filter(lambda x: self.__calibration_y.index(x) not in ind_list,self.__calibration_y)
        
    
    def delete_energies(self,ind_list):
        """
        set indexes to delete elements from energies array
        """
        if not ind_list:
            print 'no indexes to delete'
            return
        self.__energies = filter(lambda x: self.__energies.index(x) not in ind_list,self.__energies)    


    def show_points(self): #show calibration list
        if len(self.__calibration_x) == 0:
            print "Calibration peaks list: empty"
            return
        i = 0
        sort_ind = np.argsort(self.__calibration_x)
        self.__calibration_x,self.__calibration_y = np.array(self.__calibration_x)[sort_ind],np.array(self.__calibration_y)[sort_ind]
        self.__calibration_x,self.__calibration_y = list(self.__calibration_x),list(self.__calibration_y)
        print "Calibration peaks list:"
        for x,y in zip(self.__calibration_x,self.__calibration_y):
            print '( %d, x: %f, y: %f )' %(i,x,y) 
            i += 1
            
    def show_energies(self):
        if len(self.__energies) == 0:
            print "Energy list: empty"
            return
        i = 0
        print "Energy list:"
        self.__energies = sorted(self.__energies)
        for x in self.__energies:
            print '( %d, energy: %f )' %(i,x)
            i += 1


    def calibrate(self,p0 = (2,0),filename=None,ind=1):
        
        self.show_energies()
        self.show_points()         
        if len(self.__calibration_x) != len(self.__energies):
            print "Size of calibrate points array doesn't coincide the size of energy array"
            return
        
        self.__energies,self.__calibration_x = np.array(self.__energies),np.array(self.__calibration_x)
            
        def residuals(coef,y,x):
            return y - self.calibration_function(coef,self.__calibration_x)   
            
        solution = leastsq(residuals,p0, args = (self.__energies,self.__calibration_x) )
        solution = solution[0]
        fig,ax = plt.subplots()
        ax.plot(self.__calibration_x,self.__energies,'ro')
        print make_report(self.__calibration_x,solution,ax,calibration_function = self.calibration_function,filename=filename,energies=self.__energies,ind=ind)
        
        self.__energies,self.__calibration_x = list(self.__energies),list(self.__calibration_x)

                       
    def clean(self):
        lines = self.__ax.lines
        if len(lines) < 1: return
        else:
            for i in xrange(len(lines)):
                l = lines.pop(0)
                del l
        self.__calibration_x = []
        self.__calibration_y = []
        plt.show()
        
    def read(self,x,y):
        self.clean()
        try:
            self.__ax.plot(x,y,color='b',linestyle = 'steps')   
        except ValueError:
            print 'Wrong or empty arguments: x,y'
        plt.show()
                           
    def dl(self,x1, x2):
        self.__ax.set_xlim([x1,x2])
        plt.show()
        
if __name__ == '__main__':

# Simple exampe 
#    Clbr = calibrator()
#    data = np.random.randn(5000)
#    x = np.arange(-4,4,0.1)
#    y = np.histogram(data,bins=x)[0]
#    x = x[1:]
#    Clbr.read(x,y)


    from read_files import read_files, get_front_spectrs
    import os
    
    #read data
    os.chdir(r'./Files88_91')
    frame = read_files('tsn.88-tsn.91',strip_convert=True)
    hist,hsum  = get_front_spectrs(frame,tof=False,type_scale='Fchannel',visualize=True)
    hist1,hsum1 = get_front_spectrs(frame,tof=True,type_scale='Fchannel',visualize=True)
    
    #select a spectrum including alpha-peaks from events without time-of-flight mark (tof=False, first 100 channels)
    #and a peak of scattered ions which lies in area > 100 channel
    ind = 7
    def get_hist(i):
        hist_probe = hist[i]
        hist_probe[100:] = 0
        hist_probe += hist1[i]
        return hist_probe
        
    hist_probe = get_hist(ind)
    
    #create calibrator object to make calibrations by hand
    Clbr = calibrator(energies = [6040,6143,6264,6899.2,7137,7922,8699,9261], calibration_function = exponential)
        # while initiation of "calibrator" class one can set a list of energies and function for calibration of the final peak-energies list
    Clbr.read(hist_probe.index,hist_probe)
    Clbr.dl(25,500) #select an area from 25 to 500 channel
    Clbr.set_energies([7922,8699,9261,196591])
    
    print(" Select 3 most right alpha-peaks (in area <100 channel) and high peak of scattered ions (nearby 450 channel) \n \
and then put next commands into console panel to calculate calibration coefficients: \n \
    Clbr.set_energies([7922,8699,9261,196591]) #energies in MeV of the peaks \n \
    Clbr.calibrate(p0=(0,0.2,-10),filename='../Sep2014_calibrations/fission_clbr.txt',ind=ind)) #one could set initial coefficients for iterational calculating process and a filename for ouput report \n \
if you have selected some needless excess peaks, you can delete them from the list using next commands: \n \
    Clbr.show_points() # the command show a list of sorted selected peaks and their indexes \n \
    Clbr.delete_points([5,6,7]) #set the indexes of excess peaks to delete them \n \
# the same works for energies, just use commands show_energies, delete_energies in the same way.\
        ")
        #rClbr = calibrator(x,y)
        
        
        
        