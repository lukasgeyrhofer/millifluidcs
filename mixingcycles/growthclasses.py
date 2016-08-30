#!/usr/bin/env python

import numpy as np
import argparse
from scipy.stats import poisson

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

def PoissonSeedingVectors(m,n,cutoff = 1e-100,diff = False):
    if n[0] > 0:
        px = poisson.pmf(m,n[0])
        px[px<cutoff] = 0.
        px[-1] += (1. - np.sum(px))
        if diff:
            dpx = (m/n[0] - 1.)*px
    else:
        px = np.zeros(len(m))
        px[0] = 1
        if diff:
            dpx = -px
    if n[1] > 0:
        py = poisson.pmf(m,n[1])
        py[py<cutoff] = 0.
        py[-1] += (1. - np.sum(py))
        if diff:
            dpy = (m/n[1] - 1.)*py
    else:
        py = np.zeros(len(m))
        py[0] = 1
        if diff:
            dpy = -py
    if diff:
        return px,py,dpx,dpy
    else:
        return px,py



class GrowthDynamics:
    def __init__(self,growthrates = np.array([2.,1.]), yieldrates = np.array([2.,1.]), dilution = 1., mixingtime = 100., substrate = 1e4,NR_alpha = 1.,NR_precision = 1e-10, NR_maxsteps = 10000 ):
        
        self.attributes = ['growthrates','yieldrates','dilution','mixingtime','substrate','numstrains']
        self.growthrates = growthrates
        self.yieldrates = yieldrates
        self.dilution = dilution
        self.mixingtime = mixingtime
        self.substrate = substrate
        
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
      

    def checkInitialCells(self,initialcells = None):
        assert len(self.growthrates) == len(self.yieldrates)
        if not initialcells is None:
            if isinstance(initialcells,np.ndarray):
                if len(initialcells) > self.numstrains:
                    ret_ic = initialcells[:self.numstrains]
                elif len(initialcells) < self.numstrains:
                    ret_ic = np.concatenate([initialcells,np.zeros(self.numstrains - len(initialcells))])
                else:
                    ret_ic = initialcells
            else:
                ret_ic = np.ones(self.numstrains)
        else:
            ret_ic = np.ones(self.numstrains)
        return ret_ic
        

    def getGrowth(self,initialcells = None):
        ic = self.checkInitialCells(self,initialcells)
        return self.__dilution * ic * np.exp( self.__growthrates * self.__getTimeToDepletion(ic) )


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
        
    def __str__(self):
        return "growthclass object"


class StochasticGrowthDynamics(GrowthDynamics):
    
    def __init__(self,growthrates = np.array([2.,1.]), yieldrates = np.array([1.,2.]), dilution = 1., mixingtime = 100., substrate = 1e4,NR_alpha = 1.,NR_precision = 1e-10, NR_maxsteps = 10000 ):
        self.attributes = ['growthrates','yieldrates','dilution','mixingtime','substrate','numstrains']
        self.growthrates = growthrates
        self.yieldrates = yieldrates
        self.dilution = dilution
        self.mixingtime = mixingtime
        self.substrate = substrate
        
        assert len(self.growthrates) == len(self.yieldrates)        
        
        self.__NR = {'alpha':NR_alpha, 'precision2': NR_precision**2, 'maxsteps':NR_maxsteps}
        
        self.numstrains = len(self.growthrates)
        
    
    def __getNextDivision(self,population):
        totalrate = np.dot(population,self.growthrates[:len(population)])
        return np.random.choice(len(population),p = population*self.growthrates[:len(population)]/totalrate),np.random.exponential(1./totalrate)
    
    def growth(self,initialcells = None):
        n = np.array(GrowthDynamics.checkInitialCells(self,initialcells),dtype=int)
        t = 0
        s = self.substrate
        while True:
            i,dt  = self.__getNextDivision(n)
            
            if s-self.yieldrates[i] < 0: break
            if t+dt > self.mixingtime: break

            t += dt
            n[i] += 1
            s -= self.yieldrates[i]

        return min(t,self.mixingtime),n

