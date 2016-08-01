#!/usr/bin/env python

import numpy as np
import argparse


def RungeKutta4(func,xx,tt,step):
  # 4th order Runge-Kutta integration scheme
  k1 = step * func( tt        , xx )
  k2 = step * func( tt+step/2., xx+k1/2. )
  k3 = step * func( tt+step/2., xx+k2/2. )
  k4 = step * func( tt+step   , xx+k3 )
  return xx + (k1+2*k2+2*k3+k4)/6.


def addgrowthparameters(p):
    gp = p.add_argument_group(description = "Parameters for growth in droplets")
    gp.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
    gp.add_argument("-Y","--yieldrates",type=float,nargs="*",default=[1.,2.])
    gp.add_argument("-S","--substrateconcentration",type=float,default=1e4)
    gp.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
    gp.add_argument("-T","--mixingtime",type=float,default=12.)
    
    return p


class growthdynamics:
    def __init__(self,growthrates = np.array([2.,1.]), yieldrates = np.array([2.,1.]), dilution = 1., mixingtime = 100., substrate = 1e4,NR_alpha = 1.,NR_precision = 1e-10, NR_maxsteps = 10000 ):
        
        self.__attributes = ['growthrates','yieldrates','dilution','mixingtime','substrate','numstrains']
        self.__growthrates = growthrates
        self.__yieldrates = yieldrates
        self.__dilution = dilution
        self.__mixingtime = mixingtime
        self.__substrate = substrate
        
        assert len(self.__growthrates) == len(self.__yieldrates)        
        
        self.__NR = {'alpha':NR_alpha, 'precision2': NR_precision**2, 'maxsteps':NR_maxsteps}
        
    def __getTimeToDepletion(self,initialcells):
        # internal function to determine when substrate is used up
        t0 = 0
        if np.sum(initialcells) > 0.:
            # assume only single 
            t1 = max(np.log(self.__substrate*self.__yieldrates[initialcells > 0]/initialcells[initialcells > 0]+1.)/self.__growthrates[initialcells >0])
            i = 0
            while ((t1-t0)/t1)**2 > self.__NR['precision2']:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += self.__NR['alpha']*(self.__substrate-np.sum(initialcells[initialcells>0]/self.__yieldrates[initialcells>0]*(np.exp(self.__growthrates[initialcells>0]*t1)-1.)))/(np.sum(initialcells[initialcells>0]/self.__yieldrates[initialcells>0]*self.__growthrates[initialcells>0]*np.exp(self.__growthrates[initialcells>0]*t1)))
                i+=1
                # should not iterate infinitely
                if i > self.__NR['maxsteps']:
                    raise ValueError
            return min(t1,self.__mixingtime)
        else:
            return 0.
      

    def getGrowth(self,initialcells = None):
        assert len(self.__growthrates) == len(self.__yieldrates)
        if not initialcells is None:
            if isinstance(initialcells,np.ndarray):
                if len(initialcells) > self.numstrains:
                    initialcells = initialcells[:self.numstrains]
                elif len(initialcells) < self.numstrains:
                    initialcells = np.concatenate([initialcells,np.zeros(self.numstrains - len(initialcells))])
            else:
                initialcells = np.ones(self.numstrains)
        else:
            initialcells = np.ones(n)
        return self.__dilution * initialcells * np.exp( self.__growthrates * self.__getTimeToDepletion(initialcells) )


    def getGrowthMatrix(self,size):
        if isinstance(size,int):
            m = np.arange(size)
            g = np.zeros((size,size,2))
            for i in m:
                for j in m:
                    g[i,j] = self.getGrowth(initialcells = np.array([i,j]))
            return g[:,:,0],g[:,:,1]        
        elif isinstance(size,np.ndarray):
            if isinstance(size[0],np.ndarray):
                if len(size) >= 2:
                    m0 = size[0]
                    m1 = size[1]
                    g = np.zeros((len(m0),len(m1),2))
                    for i in range(len(m0)):
                        for j in range(len(m1)):
                            g[i,j] = self.getGrowth(initialcells = np.array([m0[i],m1[j]]))
                    return g[:,:,0],g[:,:,1]        
            else:
                m = size
                g = np.zeros(len(m))
                for i in range(len(m)):
                    g[i] = self.getGrowth(initialcells = np.array([m[i]]))[0]
                return g
        elif isinstance(size,(list,tuple)):
            if (len(size) >= 2) and isinstance(size[0],np.ndarray):
                m0 = size[0]
                m1 = size[1]
                g = np.zeros((len(m0),len(m1),2))
                for i in range(len(m0)):
                    for j in range(len(m1)):
                        g[i,j] = self.getGrowth(initialcells = np.array([m0[i],m1[j]]))
                return g[:,:,0],g[:,:,1]
                
    
        
    def getSingleStrainFixedPoints(self):
        t = 1./self.__growthrates * np.log(1./self.__dilution)
        y = np.array([ self.__yieldrates[i] if t[i] <= self.__mixingtime else 0. for i in range(len(self.__yieldrates))])
        return self.__dilution / (1. - self.__dilution) * self.__substrate * y
    
            
            
    def getTimeToDepletionMatrix(self,size):
        m = np.arange(size)
        t = np.zeros((size,size))
        for i in m:
            for j in m:
                t[i,j] = self.__getTimeToDepletion(initialcells = np.array([i,j]))
        return t


            
    def __getattr__(self,key):
        if key in self.__attributes:
            if   key == 'growthrates':  return self.__growthrates
            elif key == 'yieldrates':   return self.__yieldrates
            elif key == 'substrate':    return self.__substrate
            elif key == 'mixingtime':   return self.__mixingtime
            elif key == 'dilution':     return self.__dilution
            elif key == 'numstrains':   return len(self.__growthrates)

    def setGrowthRates(self,growthrates):
        if isinstance(growthrates,np.ndarray):
            self.__growthrates = growthrates
        elif isinstance(growthrates,(int,float)):
            self.__growthrates = np.array([float(growthrates)])
        elif isinstance(growthrates,(list,tuple)):
            self.__growthrates = np.array(growthrates)
    
    def setYieldRates(self,yieldrates):
        if isinstance(yieldrates,np.ndarray):
            self.__yieldrates = yieldrates
        elif isinstance(yieldrates,float):
            self.__yieldrates = np.array([yieldrates])
        elif isinstance(yieldrates,(list,tuple)):
            self.__yieldrates = np.array(yieldrates)

    def setMixingTime(self,mixingtime):
        if isinstance(mixingtime,(int,float)):
            self.__mixingtime = 1.*mixingtime
    
    def setSubstrate(self,substrate):
        if isinstance(substrate,(int,float)):
            self.__substrate = 1.*substrate
    
    def setDilution(self,dilution):
        if isinstance(dilution,(int,float)):
            self.__dilution = 1.*dilution
        


