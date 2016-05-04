#!/usr/bin/env python

import numpy as np

class growthdynamics:
    def __init__(self,growthrates = np.array([2.,1.]), yieldrates = np.array([2.,1.]), dilution = 1., mixingtime = 100., substrate = 1e4,NR_alpha = 1.,NR_precision = 1e-10, NR_maxsteps = 10000 ):
        
        self.__attributes = ['growthrates','yieldrates','dilution','mixingtime','substrate']
        self.__growthrates = growthrates
        self.__yieldrates = yieldrates
        self.__dilution = dilution
        self.__mixingtime = mixingtime
        self.__substrate = substrate
        
        assert len(self.__growthrates) == len(self.__yieldrates)        
        self.__numstrains = len(self.__growthrates)
        
        self.__NR = {'alpha':NR_alpha, 'precision2': NR_precision**2, 'maxsteps':NR_maxsteps}
        
    def __getTimeToDepletion(self,initialcells):
        # internal function to have all parameters.
        t0 = 0
        if np.sum(initialcells) > 0.:
            # initial condition for iteration is assuming only strain with highest expected yield is present
            p = (1.*initialcells/self.__yieldrates).argmax() # get ID of strain
            # set first estimate of time for single strain
            t1 = np.log(self.__substrate*self.__yieldrates[p]/initialcells[p]+1.)/self.__growthrates[p]
            i = 0
            while ((t1-t0)/t1)**2 > self.__NR['precision2']:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += self.__NR['alpha']*(self.__substrate-np.sum(initialcells/self.__yieldrates*(np.exp(self.__growthrates*t1)-1.)))/(np.sum(initialcells/self.__yieldrates*self.__growthrates*np.exp(self.__growthrates*t1)))
                i+=1
                # should not iterate infinitely
                if i > self.__NR['maxsteps']:
                    raise ValueError
            return min(t1,self.__mixingtime)
        else:
            return 0.
      

    def getGrowth(self,initialcells = None):
        assert len(self.__growthrates) == len(self.__yieldrates)        
        if initialcells is None:
            if isinstance(initialcells,np.ndarray):
                if len(initialcells) > self.__numstrains:
                    initialcells = initialcells[:self.__numstrains]
                elif len(initialcells) < self.__numstrains:
                    initialcells = np.concatenate([initialcells,np.zeros(self.__numstrains - len(initialcells))])
            else:
                initialcells = np.ones(self.__numstrains)
        else:
            initialcells = np.ones(self.__numstrains)
        return self.__dilution * initialcells * np.exp( self.__growthrates * self.__getTimeToDepletion(initialcells) )

    def getGrowthMatrix(self,size):
        m = np.arange(size)
        g = np.zeros((size,size,2))
        for i in m:
            for j in m:
                g[i,j] = self.getGrowth(initialcells = np.array([i,j]))
        return g[:,:,0],g[:,:,1]
        
        
    def getSingleStrainFixedPoints(self):
        return self.__dilution / (1. - self.__dilution) * self.__substrate * self.__yieldrates
    
    def __getattr__(self,key):
        if key in self.__attributes:
            if   key == 'growthrates':  return self.__growthrates
            elif key == 'yieldrates':   return self.__yieldrates
            elif key == 'substrate':    return self.__substrate
            elif key == 'mixingtime':   return self.__mixingtime
            elif key == 'dilution':     return self.__dilution
            
            
            
            
