#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 17:30:00 2019

Author: Max Felius

Script solving first order differential equation using the 4th order Runga-Kuttamethod.
"""
import numpy as np
import matplotlib.pyplot as plt

def ode_function(t,y1,y2):
    #variables
    m = 1 #kg
    wn = 1 #rad/s
    F0 = 1 #N
    w = 0.4*wn # rad/s
    z = 0.03
    
    #ode function to be evaluated
    feval = (F0/m)*np.sin(w*t) - 2*z*wn*y2 - wn**2 * y1
    
    return feval

def exact_solution(ts,te):
    #variables
    m = 1 #kg
    wn = 1 #rad/s
    F0 = 1 #N
    w = 0.4*wn # rad/s
    z = 0.03
    
    #initial conditions
    x0 = 0
    x_dot0 = 0    
    
    wd = wn*np.sqrt(1 - z**2)
    den = (wn**2 - w**2)**2 + (2*w*wn*z)**2
    C1 = (wn**2 - w**2)/den*F0/m
    C2 = -2*w*wn*z/den*F0/m
    A = x0*wn/wd + x_dot0/wd +(w**2 + (2*z**2 - 1)*wn**2)/den*w/wd*F0/m
    B = x0 + 2*w*wn*z/den*F0/m
    t = np.linspace(ts, te, 5000)
    x = (A*np.sin(wd*t) + B*np.cos(wd*t))*np.exp(-wn*z*t) + C1*np.sin(w*t) + C2*np.cos(w*t)
    
    plotfunction(x,t,'True Solution')
    

def RK4(ts,te,x,xdot,h):
    #RK will be evaluated at stage 4
    
    a = np.array([0.0,0.5,0.5,0.0])
    b = np.array([[0.0,0.0,0.0],[0.5,0.0,0.0],[0.0,0.5,0.0],[0.0,0.0,1.0]])
    c = np.array([1/6,1/3,1/3,1/6])
    
    t = ts
    y_inner = np.array([x,xdot])
    
    tout = [t]
    yout = [y_inner]
    
    f = np.array([[x,xdot],[0.0,0.0],[0.0,0.0],[0.0,0.0]])
    
    while t < te:
        ti = t
        yi = y_inner
        
        for i in range(4):
            t_inner = ti + a[i] * h
            y_inner = yi
            
            for j in range(i-1):
                y_inner = y_inner + h*b[i,j] * f[j,:]
            
            #ODE FUNCTION                
            y2dot = ode_function(t_inner,y_inner[0],y_inner[1])
            
            f[i,:] = np.array([y_inner[1],y2dot])
               
        h = min([h,te-t])
        t = t + h
        y_inner = yi + h*(f.T @ c)
        
        tout.append(t)
        yout.append(y_inner)
        
    return tout, yout    
    

def plotfunction(yout,tout,title):
    plt.figure()
    plt.plot(tout,yout)
    plt.grid(True)
    plt.title(title)
    
def _test():
    plt.close('all')
    
    #time vector
    ts = 0.0 #start
    te = 110 #stop (end)
    
    #initial condition
    x = 0.0
    xdot = 0.0
    
    #step size
    h = 0.01
    
    #Execute RK4 method
    tout,yout = RK4(ts,te,x,xdot,h)
    
    #plot the results
    plotfunction(np.array(yout)[:,0],tout,'RK4 Solution')
    
    exact_solution(ts,te)
    
if __name__ == '__main__':
    _test()


