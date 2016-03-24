#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

def RungeKutta4(func,xx,tt,step):
  # 4th order Runge-Kutta integration scheme
  k1 = step * func( tt       , xx )
  k2 = step * func( tt+step/2, xx+k1/2 )
  k3 = step * func( tt+step/2, xx+k2/2 )
  k4 = step * func( tt+step  , xx+k3 )
  return xx + (k1+2*k2+2*k3+k4)/6.


class dropletdynamics:
    def __init__(self,**kwargs):
        self.__x = np.array([kwargs.get('start_bacteria',10.),kwargs.get('start_nutrients',1.5),kwargs.get('start_antibiotics',0)])
        self.__antibiotics = {'zmic':kwargs.get('antibiotics_zmic',0.002),'kappa':kwargs.get('antibiotics_kappa',2.),'gamma':kwargs.get('antibiotics_gamma',10)}
        self.__nutrients = {'ks':kwargs.get('nutrients_ks',0.2),'amax':kwargs.get('nutrients_amax',0.02),'yield':kwargs.get('yieldfactor',2e5)}
        self.__time = 0
        self.__steps = 0
        self.__epsilon = kwargs.get('epsilon',0.05)
    
    def __str__(self):
        return "{:.2f} {:e} {:e} {:e}".format(self.__time,self.__x[0],self.__x[1],self.__x[2])
    
    def dxdt(self,tt,xx):
        growthrate  = self.__nutrients['amax'] * xx[1] / (self.__nutrients['ks'] + xx[1]) * (1-(1+self.__antibiotics['gamma'])*np.power(xx[2]/self.__antibiotics['zmic'],self.__antibiotics['kappa'])/(np.power(xx[2]/self.__antibiotics['zmic'],self.__antibiotics['kappa']) + self.__antibiotics['gamma']))
        yieldfactor = self.__nutrients['yield'] # still a dummy
        return np.array([growthrate * xx[0], -growthrate/yieldfactor * xx[0],0])
    
    def step(self):
        self.__x      = RungeKutta4(self.dxdt,self.__x,self.__time,self.__epsilon)
        self.__x[self.__x<0] = 0
        self.__time  += self.__epsilon
        self.__steps += 1
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-N","--start_bacteria",type=float,default=10)
    parser.add_argument("-S","--start_substrate",type=float,default=1.5)
    parser.add_argument("-B","--start_antibiotics",type=float,default=0)
    parser.add_argument("-z","--antibiotics_zmic",type=float,default=0.002)
    parser.add_argument("-k","--antibiotics_kappa",type=float,default=2.)
    parser.add_argument("-g","--antibiotics_gamma",type=float,default=10.)
    parser.add_argument("-s","--nutrients_ks",type=float,default=0.2)
    parser.add_argument("-a","--nutrients_amax",type=float,default=0.02)
    parser.add_argument("-y","--yieldfactor",type=float,default=2e5)
    parser.add_argument("-e","--epsilon",type=float,default=0.05)
    parser.add_argument("-M","--maxsteps",type=int,default=2000)
    parser.add_argument("-O","--outputstep",type=int,default=20)
    args = parser.parse_args()
    
    droplet = dropletdynamics(**vars(args))
    for i in range(args.maxsteps):
        if i%args.outputstep == 0:
            print droplet
        droplet.step()
    
        
if __name__ == "__main__":
    main()
