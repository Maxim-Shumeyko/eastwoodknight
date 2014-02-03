# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 17:44:39 2014

@author: eastwoodknight
"""

class Time:     
    
    def __init__(self,mkrsec):
              
        self.mks = mkrsec 
        self.sec = self.mks / 1000000
        self.minutes = self.sec/60
        self.hours = self.minutes/60
        self.days = self.hours/24
        
        self.mks %= 1000000
        self.sec %= 60
        self.minutes %= 60
        self.hours %= 24
        
    def __str__(self):
        return '%2d d : %2d h : %2d min : %2d sec : %6d mks' %(
                self.days,self.hours,self.minutes,self.sec,self.mks)
        
if __name__ == '__main__':
    while 1:
        a = raw_input('Enter time in mks: ')
        if not a: break
        print Time(int(a))