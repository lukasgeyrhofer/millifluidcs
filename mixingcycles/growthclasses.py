#!/usr/bin/env python

class growthdynamics:
    def __init__(self,growthrates = np.array([2.,1.]), yieldrates = np.array([2.,1.]), dilution = 1., initialcells = np.ones(2), mixingtime = np.inf, substrate = 1e4,NR_alpha = 1.,NR_precision = 1e-10, NR_maxsteps = 10000 ):
        self.__growthrates = growthrates
        self.__yieldrates = yieldrates
        self.__initialcells = initialcells
        self.__dilution = dilution
        self.__mixingtime = mixingtime
        self.__substrate = substrate
        
        self.__NR = {'alpha':NR_alpha, 'precision2': NR_precision**2, 'maxsteps':NR_maxsteps}
        
        assert len(self.__growthrates) == len(self.__yieldrates) == len(self.__initialcells)
        
    def getTimeToDepletion(self):
        # internal function to have all parameters.
        t0 = 0
        if np.sum(self.__initialcells) > 0:
            # initial condition for iteration is assuming only strain with highest expected yield is present
            p = (1.*self.__initialcells/self.__yieldrates).argmax() # get ID of strain
            # set first estimate of time for single strain
            t1 = np.log(self.__substrate*self.yieldrates[p]/self.__initialcells[p]+1.)/self.__growthrates[p]
            i = 0
            while ((t1-t0)/t1)**2 > self.__NR['precision2']:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += self.__NR['alpha']*(self.__substrate-np.sum(self.__initialcells/self.__yield*(np.exp(self.__growthrates*t1)-1.)))/(np.sum(self.__initialcells/self.__yield*self.__growth*np.exp(self.__growth*t1)))
                i+=1
                # should not iterate infinitely
                if i > self.__NR['maxsteps']:
                    raise ValueError
            return np.min(t1,self.__mixingtime)
        else:
            return 0.
      

    def getGrowth(self,initialcells = None):
        if initialcells != None:
            if isinstance(initialcells,np.ndarray):
                assert len(initialcells) == len(self.__initialcells)
                self.__initialcells = initialcells
        return self.__dilution * self.__initialcells * np.exp( self.__growthrates * self.getTimeToDepletion() )

