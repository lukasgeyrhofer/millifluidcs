#!/usr/bin/env python3

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


def AddGrowthParameters(p,allparams = False,deathrates = False,numdroplets = False,dilution = False):
    gp = p.add_argument_group(description = "Parameters for growth in droplets")
    gp.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
    gp.add_argument("-y","--yieldfactors",type=float,nargs="*",default=[1.,2.])
    if allparams or deathrates:
        gp.add_argument("-d","--deathrates",type=float,nargs="*",default=[0.,0.])
    gp.add_argument("-S","--substrateconcentration",type=float,default=1e4)
    gp.add_argument("-T","--mixingtime",type=float,default=12.)
    if allparams or dilution:
        gp.add_argument("-D","--dilution",type=float,default=2e-4)
    if allparams or numdroplets:
        gp.add_argument("-K","--numdroplets",type=int,default=1000)
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


class MicrobialStrain():
    def __init__(self,growthrate = 1.,yieldfactor = 1.,deathrate = 0.):
        self.growthrate  = growthrate
        self.yieldfactor = yieldfactor
        self.deathrate   = deathrate
    
    def __getattr__(self,key):
        if key == "growthrate":
            return self.__growthrate
        elif key == "yieldfactor":
            return self.__yieldfactor
        elif key == "deathrate":
            return self.__deathrate
    
    def __setattr__(self,key,value):
        def checkfloat(value,lowerbound = None,upperbound = None):
            try:
                checkedvalue = float(value)
            except:
                raise ValueError
            if not lowerbound is None:
                if checkedvalue < lowerbound:
                    checkedvalue = lowerbound
            if not upperbound is None:
                if checkedvalue > upperbound:
                    checkedvalue = upperbound
            return checkedvalue
        
        if key == "growthrate":
            self.__growthrate  = checkfloat(value,lowerbound = 0.)
        elif key == "yieldfactor":
            self.__yieldfactor = checkfloat(value,lowerbound = 0.)
        elif key == "deathrate":
            self.__deathrate   = checkfloat(value,lowerbound = 0.)
        else:
            super().__setattr__(key,value)

            
class Environment():
    def __init__(self,substrate = 1e4,dilution = 1.,mixingtime = 10.,numdroplets = 1000):
        self.substrate   = substrate
        self.dilution    = dilution
        self.mixingtime  = mixingtime
        self.numdroplets = numdroplets
        
    def __getattr__(self,key):
        if key == "substrate":
            return self.__substrate
        elif key == "dilution":
            return self.__dilution
        elif key == "mixingtime":
            return self.__mixingtime
        elif key == "numdroplets":
            return self.__numdroplets
    
    def __setattr__(self,key,value):
        def checkfloat(value,lowerbound = None,upperbound = None):
            try:
                checkedvalue = float(value)
            except:
                raise ValueError
            if not lowerbound is None:
                if checkedvalue < lowerbound:
                    checkedvalue = lowerbound
            if not upperbound is None:
                if checkedvalue > upperbound:
                    checkedvalue = upperbound
            return checkedvalue

        if key == "substrate":
            self.__substrate  = checkfloat(value,lowerbound = 0.)
        elif key == "dilution":
            self.__dilution   = checkfloat(value,lowerbound = 0.,upperbound = 1.)
        elif key == "mixingtime":
            self.__mixingtime = checkfloat(value,lowerbound = 0.)
        elif key == "numdroplets":
            try:
                self.__numdroplets = int(value)
            except:
                raise ValueError
            if self.__numdroplets < 1:
                self.__numdroplets = 1
        else:
            super().__setattr__(key,value)
    
    def getParams():
        return {"substrate":self.__substrate,"dilution":self.__dilution,"mixingtime":self.__mixingtime,"numdroplets":self.__numdroplets}


class GrowthDynamics:
    def __init__(self,growthrates = np.array([2.,1.]), yieldfactors = np.array([2.,1.]), deathrates = None, dilution = 1., mixingtime = 100., substrate = 1e4,NR_alpha = 1.,NR_precision = 1e-10, NR_maxsteps = 10000 ):
        
        assert len(growthrates) == len(yieldfactors)
        if deathrates is None:
            deathrates = np.zeros(len(growthrates))
        else:
            assert len(growthrates) == len(deathrates)
        
        self.strains = list()
        for a,y,d in zip(growthrates,yieldfactors,deathrates):
            self.strains.append(MicrobialStrain(growthrate = a,yieldfactor = y,deathrate = d))
            
        self.env = Environment(dilution = dilution,mixingtime = mixingtime, substrate = substrate)
        self.NR  = {'alpha':NR_alpha, 'precision2': NR_precision**2, 'maxsteps':NR_maxsteps}
    
    def addStrain(self,growthrate = 1.,yieldfactor = 1.,deathrate = 0):
        self.strains.append(MicrobialStrain(growthrate = growthrate, yieldfactor = yieldfactor, deathrate = deathrate))
    
    def delLastStrain(self):
        return self.strains.pop()
    
    def __getTimeToDepletion(self,initialcells):
        # internal function to determine when substrate is used up
        t0 = 0
        if np.sum(initialcells) > 0.:
            # assume only single strain
            t1 = max(np.log(self.env.substrate*self.yieldfactors[initialcells > 0]/initialcells[initialcells > 0]+1.)/self.growthrates[initialcells >0])
            i = 0
            while ((t1-t0)/t1)**2 > self.NR['precision2']:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += self.NR['alpha']*(self.env.substrate-np.sum(initialcells[initialcells>0]/self.yieldfactors[initialcells>0]*(np.exp(self.growthrates[initialcells>0]*t1)-1.)))/(np.sum(initialcells[initialcells>0]/self.yieldfactors[initialcells>0]*self.growthrates[initialcells>0]*np.exp(self.growthrates[initialcells>0]*t1)))
                i+=1
                # should not iterate infinitely
                if i > self.NR['maxsteps']:
                    raise ValueError
            return min(t1,self.env.mixingtime)
        else:
            return 0.
      

    def checkInitialCells(self,initialcells = None):
        try:
            ret_ic = np.array(initialcells,dtype=float)
        except:
            ret_ic = np.ones(self.numstrains,dype=float)
        if len(ret_ic) < self.numstrains:
            ret_ic = np.concatenate((ret_ic,np.zeros(self.numstrains - len(ret_ic),dtype=float)))
        elif len(ret_ic) > self.numstrains:
            ret_ic = ret_ic[:self.numstrains]
        return ret_ic
        

    def Growth(self,initialcells = None):
        ic = self.checkInitialCells(initialcells)
        return self.env.dilution * ic * np.exp( self.growthrates * self.__getTimeToDepletion(ic) )


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
        t = 1./self.growthrates * np.log(1./self.env.dilution)
        y = np.array([ self.yieldfactors[i] if t[i] <= self.mixingtime else 0. for i in range(self.numstrains)])
        return self.env.dilution / (1. - self.env.dilution) * self.env.substrate * y
    
            
            
    def getTimeToDepletionMatrix(self,size):
        m = np.arange(size)
        t = np.zeros((size,size))
        for i in m:
            for j in m:
                t[i,j] = self.__getTimeToDepletion(initialcells = np.array([i,j]))
        return t

    
    def __getattr__(self,key):
        if key == "numstrains":
            return len(self.strains)
        elif key == "growthrates":
            return np.array([self.strains[i].growthrate for i in range(self.numstrains)])
        elif key == "yieldfactors":
            return np.array([self.strains[i].yieldfactor for i in range(self.numstrains)])
        elif key == "deathrates":
            return np.array([self.strains[i].deathrate for i in range(self.numstrains)])
        

    def __setattr__(self,key,value):
        if key == "growthrates":
            try:
                tmp = np.array(value,dtype=float)
            except:
                raise ValueError
            assert len(tmp) == self.numstrains
            for i in range(self.numstrains):
                self.strains[i].growthrate = tmp[i]
        elif key == "yieldfactors":
            try:
                tmp = np.array(value,dtype=float)
            except:
                raise ValueError
            assert len(tmp) == self.numstrains
            for i in range(self.numstrains):
                self.strains[i].yieldfactor = tmp[i]
        elif key == "deathrates":
            try:
                tmp = np.array(value,dtype=float)
            except:
                raise ValueError
            assert len(tmp) == self.numstrains
            for i in range(self.numstrains):
                self.strains[i].deathrate = tmp[i]
        else:
            super().__setattr__(key,value)
        
    def setMixingTime(self,mixingtime):
            self.env.mixingtime = mixingtime
    def setSubstrate(self,substrate):
            self.env.substrate = substrate
    def setDilution(self,dilution):
            self.env.dilution = dilution


class StochasticGrowthDynamics(GrowthDynamics):
    def __init__(self,**kwargs):
        GrowthDynamics.__init__(self,kwargs)
        self.__lastgrowthtime = np.nan

    def __getNextDivision(self,population):
        totalrate = np.dot(population,self.growthrates[:len(population)])
        return np.random.choice(len(population),p = population*self.growthrates[:len(population)]/totalrate),np.random.exponential(1./totalrate)
    
    def checkInitialCells(self, initialcells = None):
        return np.array(GrowthDynamics.checkInitialCells(initialcells),dtype=int)
    
    def Growth(self,initialcells = None):
        n = self.checkInitialCells(self,initialcells)
        t = 0
        s = self.env.substrate
        while True:
            i,dt  = self.__getNextDivision(n)
            
            if s-self.yieldfactor[i] < 0: break
            if t+dt > self.mixingtime: break

            t += dt
            n[i] += 1
            s -= self.yieldrates[i]
        self.__lastgrowthtime = min(t,self.env.mixingtime)
        return n
    
    def __getattr__(self,key):
        if key == "lastgrowthtime":
            if self.__lastgrowthtime is np.nan:
                raise ValueError("StochasticGrowthDynamics.Growth(initialcells) was not yet called")
            else:
                return self.__lastgrowthtime
        else:
            super().__getattr__(self,key)
