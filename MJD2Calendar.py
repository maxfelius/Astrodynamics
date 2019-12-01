#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Max Felius
studentnumber: 4250788
email: maxfelius@hotmail.com

AE4878 - MGOD

Assignment 2 - BASICS 3

TODO:
    -input in class should either be F or HH:MM:SS
    -Output seconds should have more decimal points
    -make the classes also accept lists and tuples?
"""

#imports
import numpy as np

#constants

#Modified Julian Day Object
class MJD(object):
    def __init__(self,mjd):
        self.mjd = mjd
    
    def __str__(self):
        return 'Modified Julian day is: {}'.format(self.mjd)

    def printCalendarDate(self):
        obj = self.MJD2CalDate()
        print('Converted Calendar Date: Y={}, M={}, D={}, HH={}, MM={}, SS={}'.format(obj[0],obj[1],obj[2],obj[3],obj[4],obj[5]))
    
    def F2hhmmss(self,F):
        time1 = 24*F
        time2 = 24*F-int(24*F)
        
        hh = np.fix(time1)

        mm = np.fix(time2*60) #minutes
        
        ss = np.fix((time2*60-mm)*60) #seconds
        return hh,mm,ss
        
    def MJD2CalDate(self):
        JD = self.mjd + 2400000.5
        JD0 = JD+0.5
        
        F = JD0-int(JD0)
        
        L = np.fix(JD0 + 68569)
        N = np.fix(4*L/146097)
        
        L = L - np.fix((146097*N + 3)/4)
        
        I = np.fix((4000*(L+1))/1461001)
        
        L = L - np.fix(1461 * (I/4)) + 31
        
        J = np.fix((80*L)/2447)
        
        #day
        D = L - np.fix(2447* (J/80))
        
        L = np.fix(J/11)
        
        #month
        M = J + 2 - 12*L
        #year
        Y = 100*(N-49) + I + L
        
        #time conversion
        hh,mm,ss = self.F2hhmmss(F)
        
        return Y,M,D,hh,mm,ss

class CalendarDate(object):
    def __init__(self,Y,M,D,F):
        self.Y = Y
        self.M = M
        self.D = D
        self.F = F
        
    def __str__(self):
        hh,mm,ss = self.F2hhmmss(self.F)
        return 'Calender Date is: Y={}, M={}, D={}, HH={}, MM={}, SS={}'.format(self.Y,self.M,self.D,hh,mm,ss)

    def printMJD(self):
        print('Converted Modified Julian Day: {}'.format(self.CalDate2MJD()))
            
    def F2hhmmss(self,F):
        time1 = 24*F
        time2 = 24*F-int(24*F)
        
        hh = np.fix(time1)

        mm = np.fix(time2*60) #minutes
        
        ss = np.fix((time2*60-mm)*60) #seconds
        return hh,mm,ss

    def CalDate2MJD(self):
        C = np.fix((self.M-14)/12)
        
        JD0 = self.D - 32075 + np.fix(1461*(self.Y + 4800 + C)/4) \
        + np.fix(367.0*(self.M - 2 - C * 12)/12) \
        - np.fix(3*(np.fix((self.Y + 4900 + C))/100)/4)
        
        JD = JD0 + self.F -0.5
        
        
        MJD = JD - 2400000.5
        
        return MJD
    
def hhmmss2F(hh,mm,ss):
    F = (hh/24) + (mm/(24*60)) + (ss/(24*3600))
    return F

def _BASICS3():
    MJD_list = [49730.238000,51601.498000,51604.600000,53770.998000]
    
    CalDate = [(2000,1,18,13,15,16.668),(2000,2,28,9,9,23.800),
               (2000,2,29,18,59,1.001),(2008,10,4,1,1,1.001)]
    
    for item in MJD_list:
        print(MJD(item))
        MJD(item).printCalendarDate() 
        
    for item in CalDate:
        F = hhmmss2F(item[3],item[4],item[5])        
        print(CalendarDate(item[0],item[1],item[2],F))
        CalendarDate(item[0],item[1],item[2],F).printMJD() 
    

def _test():
    MJD1 = 54132.25
    
    print(MJD(MJD1))
    MJD(MJD1).printCalendarDate() #Y=2007, M=2, D=1, 


    # calendar to MJD
    Y1 = 2007
    M1 = 2
    D1 = 1
    F1 = 0.25
    
    print(CalendarDate(Y1,M1,D1,F1))
    CalendarDate(Y1,M1,D1,F1).printMJD() #54132.25

if __name__ == '__main__':
    _test()
    print('\nBASICS 3 assingment 2 - Results...\n')
    _BASICS3()