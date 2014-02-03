# -*- coding: utf-8 -*-
"""
###################################################################
class Visualise
__init__(self,data) # data - list, np.array, etc. 1Dimension
__call__ # show calibration menu
"""
import matplotlib.pyplot as plt
import numpy as np
#from matplotlib.widgets import Button
from hist_probe import fit_commander
from scipy.optimize import leastsq
from scipy import linalg
from scipy import poly1d
#from scipy.optimize import minimize
#import matplotlib.axes as mpax

class Visualise:
    """
    ##############################
    Allow to make calibrations and fitting with data histograms
    in a simple visualize form. Put in your data, 
    simple click twice LMB to choose an area and then click RMB to calibrate.
    ##############################
    1-st step: initialize
    __init__(self,data)
    data <- ndarray[X,Y]
    example: clbr = Visualize(spectr)
    2-st: start work window
    example: clbr()
    3-st: use the array with calibration results in .buf
    clbr.buf <- array of {[parameters]: [values]}
    # fit_result = dict(zip(parameters,solution[0]))
    examle: 
        for dicts in clbr.buf:
        print '---------------'
        for key in dicts:
            print '%5s:  %2.2f' % (key, dicts[key] )
    
    """
    def __init__(self, data, show = True):
        # ñîçäàþ îêíà äëÿ ãðàôèêîâ
        self.buf = [] # буфер для накопления отфитированных точек
        if show:
            # cоздание окон и холста
            self.fig = plt.figure()
            self.ax1 = plt.subplot(211)
            #self.ax2 = plt.subplot(232)
            self.ax3 = plt.subplot(212)
            self.ax1.legend()
            # íîìåð ãðàôèêà èç ñïèñêà, êîòîðûé áóäåò îòîáðàæàòüñÿ ïðè çàïóñêå
            # создание кнопок и графиков, назначение функций для кнопок и инициализация переменных
            self.data = data
        
            self.series3, = self.ax3.plot( np.array(data[0]),np.array(data[1]), linestyle = 'steps')
            # íàçàíà÷åíèå ïðîöåäóðû äëÿ îáðàáîòêè äàííûõ â îêíå ñïåêòðîâ
            # назначение процедуры выделения и фитирования области 
            br = fit_commander(self.ax3,self.data[1],
                           self.data[0],#"""np.arange(1,len(self.data[self.strip])+1)""",
                           self.fit )
            self.fig.canvas.mpl_connect('button_press_event', lambda event:br.onclick(event) ) 
                                    #lambda event:fit_commander.onclick(br,event))#br.onclick)
        
        
    def __call__(self):
        plt.show()
        
    """
    def redraw(self):
        self.ax3.cla()
        self.ax3.set_title( str(self.strip+1)+' / '+str(self.strip_count) )
        self.ax3.plot( np.arange(1,len(self.data[0]) + 1),self.data[self.strip], linestyle = 'steps')
        plt.draw()
        br = fit_commander(self.ax3,self.data[self.strip],
                           np.arange(1,len(self.data[self.strip])+1), self.fit_method )
        self.fig.canvas.mpl_connect('button_press_event', lambda event:br.onclick(event) ) 
                                    #lambda event:fit_commander.onclick(br,event))#br.onclick)
       
    def next_click(self,event):
        self.strip += 1
        if self.strip == self.strip_count:
            self.strip = 0
        self.redraw()

    def prev_click(self,event):
        self.strip -= 1
        if self.strip < 0:
            self.strip = self.strip_count-1
        self.redraw()

    def fit_method(self,x,y): 
        
        def func(p,x):
            return np.exp( -((x - p[0])/p[1])**2 /2 )/np.sqrt( 1 / np.pi) / p[1]

        def residuals(p,y,x):
            return y - func(p,x)

        def Jac(p0,y,x):
            dfdx0 = -( x - p0[0] )*func(p0,x)/2/p0[1]**2 #* residuals(p0,y,x) 
            dfdsigm = (1/p0[1] - ( x - p0[0])**2 /(p0[1]**3) /2 )*func(p0,x)#* residuals(p0,y,x) 
            deriv = np.c_[dfdx0, dfdsigm]
            return deriv

        p0 = ( float(np.sum(x*y))/np.sum(y),np.std(x))
        print p0
        #âû÷èñëåíèå ïàðàìåòðîâ ôèòà ìåòîäîì Ëåâåíáåðãà-Ìàðêâàðäòà
        solution = leastsq(residuals, p0, Dfun = Jac,
                   args = (y,x) )
        self.ax1.plot(x,y,'b', linestyle = 'steps',label = 'raw data')
        self.ax1.plot(x, func(solution[0],x),'r--',label = 'fit graph')
        print solution[0]
        return solution[0] # âîçâðàùàåò ïàðàìåòðû ôèòà 
    """
    def fit_method(self,x,y):
        p0 = ( float(np.sum(x*y))/np.sum(y),np.std(x),np.sum(y)/np.sqrt(np.std(x))/np.sqrt(1/(2*np.pi)),
               float(y[-1] - y[0])/(x[-1] - x[0]), 1 )
        
        def func(p,x):
            # P: 0 - x0; 1 - sigm; 2 - A; 3 - B; 4 - C
            return ( p[2]* np.exp( -((x - p[0])/p[1])**2 /2 ) /np.sqrt( 1/(2*np.pi))/ p[1] 
                   +p[3]*(p[0]*1.003 - x) + p[4] )
        
        def residuals(p,y,x):
            return y - func(p,x)

        #derivative
        def Jac(p0,y,x):
            dfdx0 = -( x - p0[0] )*func(p0,x)/2/p0[1]**2 * residuals(p0,y,x) 
            dfdsigm = (1/p0[1] -( x - p0[0])**2 /(p0[1]**3) /2 )*func(p0,x)
            dfdA = - np.exp( -((x - p0[0])/p0[1])**2 /2 ) /np.sqrt( 1/(2*np.pi))/ p0[1]
            dfdB = - p0[0]*1.003 + x
            dfdC = - np.ones(len(x),dtype = float)
            deriv = np.c_[dfdx0, dfdsigm,dfdA,dfdB,dfdC]
            return deriv

        solution = leastsq(residuals, p0, Dfun = Jac,
                   args = (y,x) )
        return solution
        
        
    def fit(self,x,y):

        solution = self.fit_method(x,y)
        self.ax1.cla()
        self.ax1.plot(x,y,'b', linestyle = 'steps',label = 'raw data')
        self.ax1.plot(x, func(solution[0],x),'r--',label = 'fit graph')
        plt.legend()
        plt.draw()
        # P: 0 - x0; 1 - sigm; 2 - A; 3 - B; 4 - C
        parameters = ['x0','sigm','A','B','C']
        fit_result = dict(zip(parameters,solution[0]))
        self.buf.append( fit_result )
         
      
    def calibrate(self,msv_clbr = [4824,5155,5499]):
        if self.buf == []:
            return 'empty buffer. Try to fit some peaks.'
        else:
            xi = []
            for i in self.buf:
                xi.append( i['x0'] )
            xi = np.array(xi)
            A = np.c_[np.ones(len(xi)),xi,xi**2]
            solution = linalg.lstsq(A,msv_clbr)
            solution = solution[0].tolist()
            solution.reverse()
            pol = poly1d(solution)
            plt.plot(xi,msv_clbr,'x',
                     xi,pol(xi),'r')
            return solution
            
            
    def autocalibrate(self):
        h = float(self.data[1].max())/2
        msv = np.array(self.data[1]) > h
        flag = False
        clbr_series = []
        buf = []
        j = 0
        
        for i in msv:
            #verify if inner subcequent is finished
            if not i and flag:
                flag = False
                if len(buf) >= 4:
                    clbr_series.append(buf)
                    buf = []
                else:
                    buf = []
            flag = i
            if i:
                buf.append( [self.data[0][j],self.data[1][j]] )   
            j += 1
        
        if len(buf) >= 4:
                    clbr_series.append(buf)
        
        #gluing "holes" in the peaks array
        series = [] 
        flag = False        
        for i in range(1,len(clbr_series)):
            if clbr_series[i-1][-1][0] + 10 >= clbr_series[i][0][0]:
                if flag:
                    series.pop()
                    flag = False
                series.append(clbr_series[i-1]+clbr_series[i])
            else:
                if i == 1:
                    series.append(clbr_series[0])
                series.append(clbr_series[i])
                flag = True
        #find +/- 4sigm intervals around the peaks and send them for fitting
        indexes = []
        solutions = []
        for i in series:
            x = np.array(i)[:,0]
            y = np.array(i)[:,1]
            print x,'\n',y
            mean = (x*y).sum()/y.sum()
            sigm = len(x)#x.std()
            print 'mean,sigm: ',mean,sigm
            index = (self.data[0] > mean - 5*sigm )&(self.data[0] < mean + 5*sigm)
            indexes.append(index)
            X = self.data[0][index]
            Y = self.data[1][index]
            print X,'\n',Y
            solution = self.fit_method(X,Y)
            solutions.append(solution[0])
            
        def func(p,x):
            # P: 0 - x0; 1 - sigm; 2 - A; 3 - B; 4 - C
            return ( p[2]* np.exp( -((x - p[0])/p[1])**2 /2 ) /np.sqrt( 1/(2*np.pi))/ p[1] 
                   +p[3]*(p[0]*1.003 - x) + p[4] )
        
        for i,j in zip(indexes,solutions):
            x = self.data[0][i]#np.array(i)[:,0]
            y = func(j,x)
            self.ax3.plot(x,y,'-r')
        plt.show()
    
        self.buf += solutions
        return solutions
                                    
        
         
          
        

        
if __name__ == '__main__':
    spectr = []
    #x = np.random.randn(100000);
    def func(p,x):
        # P: 0 - x0; 1 - sigm; 2 - A; 3 - B; 4 - C
        return  p[2]* np.exp( -((x - p[0])/p[1])**2 /2 ) /np.sqrt( 1/(2*np.pi))/ p[1]

    p = (20,5,800,0,0)
    x = np.linspace(0,40)
    y = func(p,x)
    spectr.append(x)
    spectr.append( y )#np.histogram(x,20)[0] );
    #spectr.append( np.sin(np.arange(500))+1 )
    #print spectr
    a = Visualise(spectr)
    a()
    print '\nFit results: '
    for dicts in a.buf:
        print '---------------'
        for key in dicts:
            print '%5s:  %2.2f' % (key, dicts[key] )


        
        
