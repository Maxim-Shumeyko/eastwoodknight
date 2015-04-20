# -*- coding: cp1251 -*-
''' Модуль предназначен для выделения области гистограммы и последующей
передачи этих данных для обработки. You can use this module to dole out part of
a histogramm to pass the data for processing. 

Пример использования:
    from mpl_cmd import fit_commander
    fig,ax = plt.subplots()
    data = np.random.randn(800)  
    y,x,patches = ax.hist(data,30,picker = 5) 
    # в качестве аргумента передаются (self,ax,y,x, method,method2=None):
        система координат (axes)
        данные о гистограмме y, x
        метод обработки данных ( method(datax,datay) ) для ЛКМ
        доп метод обработки для средней кнопки СКМ
    # формат данных (datax,datay), type(datax) == []
    def method(x,y): print x[1:],y 
    br = fit_commander(ax,y,x,method)
    fig.canvas.mpl_connect('button_press_event',
                            lambda event: evbr.onclick(event) )
    plt.show()
    # выделение области левой кнопкой мыши (press LMB - to dole out)
    # передача на обработку - правой(press RMB - for processing)
''' 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class fit_commander:
    
    def __init__(self,ax,method,method2=None):
        self.ax = ax
        self.clicks = 0
        self.span = self.ax.axvspan(0,0, facecolor = 'b', alpha = 0.2, visible = True)
        self.datax, self.datay = [],[]
        self.method = method
        self.method2 = method2
        
    def onclick(self, event):
        if len(self.ax.lines) == 0:
            return
        else:
            line = self.ax.lines[0]
            self.x = line.get_xdata()
            self.y = line.get_ydata(True)
        # обработать все данные гистограммы
        if event.dblclick:
            if self.method2:
                    self.method2(self.x,self.y)
                    return
            self.method(self.x,self.y )
            return
        
        # выделить область
        if event.button == 1:
            
            # левая граница
            if self.clicks == 0:
                self.x1 = event.xdata
               
                
            # правая граница и выделение области
            elif self.clicks == 1:                
                self.x2  = event.xdata 
                
                #fig.canvas.draw()

                def f_ind(x,x1):
                    if x1 < x[1]:
                        return 0#self.y[0]
                    elif x1 > x[-1] :
                        return len(x) + 1
                    for i in x:
                        if x1 < i: break
                    return x.tolist().index(i)- 1
                
                if self.x1 > self.x2:
                    self.x1, self.x2 = self.x2, self.x1

                if self.x1 >= self.x[0] and self.x2 >= self.x[0]:
                    i1 = f_ind(self.x,self.x1)
                    i2 = f_ind(self.x,self.x2) 
                    if self.x[i2] < self.x[-1]:
                        i2 += 1
                else:
                    self.clicks -= 1
                    return
                self.x1 = self.x[i1]
                self.x2 = self.x[i2]
                self.span.set_xy( [ [self.x1,0], [self.x1,1], [self.x2,1],
                                [self.x2,0], [self.x1,0] ] )
                self.datay = self.y[i1:i2]
                self.datax = self.x[i1:i2]
                
            self.clicks = (self.clicks + 1)%2

        # обработка данных в выделенной области по правому клику        
        if event.button == 3:
            self.span.set_xy( [ [0,0], [0,1], [0,1],
                                [0,0], [0,0] ] )
            #fig.canvas.draw()
            self.method(self.datax,self.datay )
        plt.draw()
        
if __name__ == '__main__':
    fig,ax = plt.subplots()
    data = np.random.randn(800)  
    y = np.histogram(data,bins = np.linspace(-4,4,20) )[0]
    x = np.linspace(-4,4,19)
    ax.plot(x,y,linestyle = 'steps')
    def method(x,y): print x,y        
    br = fit_commander(ax,method)
    fig.canvas.mpl_connect('button_press_event',
                           lambda event:br.onclick(event))
    plt.show()
