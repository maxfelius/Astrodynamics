#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 19:07:04 2019

Author: Max Felius
"""

"""Imported libraries"""
import numpy as np
import math

"""Global Constants"""
mu = 398600.441*10**9 #m3/s2 - Gravitational parameter Earth

class VectorOperators(object):
    ##https://en.wikipedia.org/wiki/Euclidean_vector#Basic_properties
    def __init__(self):
        """Constructor           Initiate self
        """
        
    def NormVector(self,v):
        """Input is a 1x3 vector"""
        assert len(v) == 3
        return np.sqrt(v[0]**2+v[1]**2+v[2]**2)
    
    def Normalisation(self,v):
        """input 1x3 vector"""
        return v/self.NormVector(v)   

class CartesianObject(object):
    def __init__(self,x,y,z,xdot,ydot,zdot):
        """constructor"""
        self.x = x
        self.y = y
        self.z = z
        self.xdot = xdot
        self.ydot = ydot
        self.zdot = zdot
        
    def __str__(self):
        """Returns Cartesian coordinates of the object"""
        return 'Cartesian Point: (x,y,z,xdot,ydot,zdot) = ({},{},{},{},{},{})'.format(self.x,self.y,self.z,self.xdot,self.ydot,self.zdot)
    
    def __hash__(self):
        """The hash"""
        return hash((self.x,self.y,self.z,self.xdot,self.ydot,self.zdot))
    
    def PrintKepler(self,angletype):
        print('Kepler elements: (a,e,i,RAAN,omega,Theta) = {}'.format(self.Cart2Kep(angletype)))
    
    def semimajor(self,r,V):
        return 1/((2/r)-(V**2/mu))
    
    def eccentricity(self,V,h,r):
        
        e_temp = ((np.cross(V,h))/mu)-(r/VectorOperators().NormVector(r))

        #checking if the orbit is an circle or an ellipse
        assert VectorOperators().NormVector(e_temp) < 1
        
        return e_temp
        
    def inclination(self,h):
        return math.acos(h[2]/VectorOperators().NormVector(h))
    
    def RAAN(self,N):
        
        Nx = N[0]
        Ny = N[1]
        
        N_vec = [Nx,Ny,0]
        
        Nxy = VectorOperators().NormVector(N_vec)
        
        return math.atan2(Ny/Nxy,Nx/Nxy)
    
    def omega(self,e,r,N,h):
        Nhat = VectorOperators().Normalisation(N) #N/VectorOperators.NormVector(N)
        
        if  np.cross(Nhat,e)@h.T > 0:
            sign = 1
        else:
            sign = -1
            
        ehat = VectorOperators().Normalisation(e) #e/VectorOperators.NormVector(e)    
            
        return sign * math.acos(ehat@Nhat.T)
    
    def theta(self,e,r,h):
        if (np.cross(e,r)@h.T) > 0:
            sign = 1
        else:
            sign = -1
        
        ehat = VectorOperators().Normalisation(e) #e/VectorOperators.NormVector(e)   
        rhat = VectorOperators().Normalisation(r) #r/VectorOperators.NormVector(r)
        
        return sign*math.acos(rhat@ehat.T)
    
    #Return individual values
    def a(self,angletype):
        print(self.Cart2Kep(angletype)[0])
    def e(self,angletype):
        print(self.Cart2Kep(angletype)[1])
    def i(self,angletype):
        print(self.Cart2Kep(angletype)[2])
    def AscendingNode(self,angletype):
        print(self.Cart2Kep(angletype)[3])
    def perigee(self,angletype):
        print(self.Cart2Kep(angletype)[4])
    def TrueAnomaly(self,angletype):
        print(self.Cart2Kep(angletype)[5])    
        
    def Cart2Kep(self,angletype):   
        #create radius and velocity vector
        r = (self.x,self.y,self.z);
        V = (self.xdot,self.ydot,self.zdot)
        
        #calculate the norm for the radius and velocity vector
        r_norm = VectorOperators().NormVector(r)
        V_norm = VectorOperators().NormVector(V)
        
        #Calculate angular momentum
        h = np.cross(r,V) #correct
        
        #calculate Kepler elements
        a_out = self.semimajor(r_norm,V_norm)
        e_out = self.eccentricity(V,h,r)
        i_out = self.inclination(h)
        
        N = [-h[1],h[0],0]
        
        OMEGA_out = self.RAAN(N)
        
        om_out = self.omega(e_out,r,N,h)
        
        tet_out = self.theta(e_out,r,h)
        
        #Check if the output should be in radians or degrees
        if angletype == 'rad':
            return a_out,VectorOperators().NormVector(e_out),i_out,OMEGA_out,om_out,tet_out
        elif angletype == 'deg':
            return a_out,VectorOperators().NormVector(e_out),i_out*180/math.pi,OMEGA_out*180/math.pi,om_out*180/math.pi,tet_out*180/math.pi
        else:
            print('Please specify an angletype.')
    
    """ Not enough time to this part of the script and not necessary
    
    def MeanAnomaly(self,r,a,e):
        #a is a scaler, r and e are vectors
        E = self.EccentricityAnomaly(r,a,e)
        return E - VectorOperators().NormVector(e) *math.sin(E)
        
    def EccentricityAnomaly(self,r,a,e):
        #a is a scaler, r and e are vectors
        return math.acos((VectorOperators().NormVector(r)/(-VectorOperators().NormVector(e)*a))-1)
    
    #print anomalies
    def E(self,r,a,e):
        #a is a scaler, r and e are vectors
        print(self.EccentricityAnomaly(r,a,e))
    
    def M(self,r,a,e):
        print(self.MeanAnomaly(r,a,e))
        
    """
        
    
class KeplerObject(object):
    def __init__(self,a,e,i,RAAN,omega,theta):
        """constructor"""
        self.a = a
        self.e = e
        self.i = i 
        self.RAAN = RAAN
        self.omega = omega
        self.theta = theta
        
        #check if it is an circular of elliptical orbit
        assert self.e < 1

    def __str__(self):
        """Returns Kepler elements of the object"""
        return 'Kepler elements: (a,e,i,RAAN,omega,theta) = ({},{},{},{},{},{})'.format(self.a,self.e,self.i,self.RAAN,self.omega,self.theta)
    
    def __hash__(self):
        """The hash"""
        return hash((self.a,self.e,self.i,self.RAAN,self.omega,self.M))
    
    def PrintCartesian(self):
        print('Cartesian Coordinates: (x,y,z,xdot,ydot,zdot) = {}'.format(self.Kep2Cart()))
    
    #Return individual values
    def x(self):
        print(self.Kep2Cart()[0])
    def y(self):
        print(self.Kep2Cart()[1])
    def z(self):
        print(self.Kep2Cart()[2])
    def xdot(self):
        print(self.Kep2Cart()[3])
    def ydot(self):
        print(self.Kep2Cart()[4])
    def zdot(self):
        print(self.Kep2Cart()[5])
    
    def Kep2Cart(self):
        l1 = math.cos(self.RAAN)*math.cos(self.omega) - \
            math.sin(self.RAAN)*math.sin(self.omega)*math.cos(self.i)
        l2 = -math.cos(self.RAAN)*math.sin(self.omega) - \
            math.sin(self.RAAN)*math.cos(self.omega)*math.cos(self.i)
        
        m1 = math.sin(self.RAAN)*math.cos(self.omega) + \
            math.cos(self.RAAN)*math.sin(self.omega)*math.cos(self.i)
        m2 = -math.sin(self.RAAN)*math.sin(self.omega) + \
            math.cos(self.RAAN)*math.cos(self.omega)*math.cos(self.i)
        
        n1 = math.sin(self.omega)*math.sin(self.i)
        n2 = math.cos(self.omega)*math.sin(self.i)
        
        #for an ellipse
        E = 2*math.atan(math.tan(self.theta/2)*math.sqrt((1-self.e)/(1+self.e)))
        
        r = self.a*(1-self.e*math.cos(E))
        
        ski = r*math.cos(self.theta)
        nu = r*math.sin(self.theta)
        
        #Compute x,y,z coordinates
        x_out = l1*ski + l2*nu
        y_out = m1*ski + m2*nu
        z_out = n1*ski + n2*nu
        
        #
        H = math.sqrt(mu*self.a*(1-self.e**2))

        #Compute xdot, ydot and zdot
        xdot_out = (mu/H)*(-l1*math.sin(self.theta) + \
                l2*(self.e + math.cos(self.theta)))
            
        ydot_out = (mu/H)*(-m1*math.sin(self.theta) + \
                m2*(self.e + math.cos(self.theta)))
        
        zdot_out = (mu/H)*(-n1*math.sin(self.theta) + \
                n2*(self.e + math.cos(self.theta)))
        
        return x_out,y_out,z_out,xdot_out,ydot_out,zdot_out
        
def _test():
    print('\nTest 1\n')
    x = -2700816.14 #m
    y = -3314092.80 #m
    z = 5266346.42 #m
    xdot = 5168.606550 #m/s
    ydot = -5597.546618 #m/s
    zdot = -868.878445 #m/s
    
    CartesianObject(x,y,z,xdot,ydot,zdot).PrintKepler('deg')

    a = 6787746.891
    e = 0.000731104
    i = 51.68714486*(math.pi/180) # 0.9021108587697676
    RAAN = 127.5486706*(math.pi/180) # 2.2261442580764985
    omega = 74.21987137*(math.pi/180) # 1.2953811258454424
    theta = 24.10027677*(math.pi/180) # 0.4206291802671584
    
    KeplerObject(a,e,i,RAAN,omega,theta).PrintCartesian()
    
def _test2():
    print('\nTest 2\n')
    x = 3126974.99
    y = -6374445.74
    z = 28673.59
    xdot = -254.91197
    ydot = -83.30107
    zdot = 7485.70674
    
    CartesianObject(x,y,z,xdot,ydot,zdot).PrintKepler('deg')
    
    a = 7096137.00
    e = 0.0011219
    i = 92.0316*(math.pi/180)
    RAAN = (296.1384)*(math.pi/180)
    omega = 120.6878*(math.pi/180)
    theta = (239.5437)*(math.pi/180)
    
    KeplerObject(a,e,i,RAAN,omega,theta).PrintCartesian()

def _Basics1():
    print('\nBasics-1\n')
    #assignment 1
    x = 8751268.4691 #m
    y = -7041314.6869 #m
    z = 4846546.9938 #m
    xdot = 332.2601039 #m/s
    ydot = -2977.0815768 #m/s
    zdot = -4869.8462227 #m/s
    
    CartesianObject(x,y,z,xdot,ydot,zdot).PrintKepler('deg')
    
    #assignment 2
    a = 12158817.9615 #m
    e = 0.014074320051 # [-]
    i = 52.666016957*(math.pi/180) # radians
    RAAN = 323.089150643*(math.pi/180) # radians
    omega = 148.382589129*(math.pi/180) # radians
    M = 112.192638384*(math.pi/180) # radians

    KeplerObject(a,e,i,RAAN,omega,M).PrintCartesian()

if __name__ == '__main__':
    _test()
    _test2()
    _Basics1()