#!/usr/bin/env python3
#-*- coding: utf-8 -*-


'''
==================================
=  growthclasses.py
==================================

  Contains code that generates the dynamics during a
  single cycle of the droplet dilution. Both, deterministic
  and stochastic, growth types are included as their own
  classes.
  
  In addition, several helper routines help with heavily reused
  code, to reduce the length of actual scripts that compute
  interesting properties of the dynamics.


  Lukas Geyrhofer, l.geyrhofer@technion.ac.il, 2016

'''



import numpy as np
import argparse
from scipy.stats import poisson
import itertools

def RungeKutta4(func,xx,tt,step):
  # 4th order Runge-Kutta integration scheme
  k1 = step * func( tt        , xx )
  k2 = step * func( tt+step/2., xx+k1/2. )
  k3 = step * func( tt+step/2., xx+k2/2. )
  k4 = step * func( tt+step   , xx+k3 )
  return xx + (k1+2*k2+2*k3+k4)/6.


def AddGrowthParameters(p,allparams = False,deathrates = False,numdroplets = False,dilution = False,
                        defaultgrowthrates = [2.,1.],defaultyieldfactors = [1.,2.],defaultdeathrates = [0.,0.],
                        defaultsubstrate = 1e4, defaultmixingtime = None,defaultdilution = 2e-4, defaultnumdroplets = 1000):
    # Helper routine to generate all cmdline parameters for microbial growth
    gp = p.add_argument_group(description = "Parameters for growth in droplets")
    gp.add_argument("-a","--growthrates",type=float,nargs="*",default=defaultgrowthrates)
    gp.add_argument("-y","--yieldfactors",type=float,nargs="*",default=defaultyieldfactors)
    if allparams or deathrates:
        gp.add_argument("-d","--deathrates",type=float,nargs="*",default=defaultdeathrates)
    gp.add_argument("-S","--substrateconcentration",type=float,default=defaultsubstrate)
    gp.add_argument("-T","--mixingtime",type=float,default=defaultmixingtime)
    if allparams or dilution:
        gp.add_argument("-D","--dilution",type=float,default=defaultdilution)
    if allparams or numdroplets:
        gp.add_argument("-K","--numdroplets",type=int,default=defaultnumdroplets)
    return p


def AddNRParameters(p):
    # Helper routine to generate cmdline parameters to change default behaviour of NR iterations
    nrp = p.add_argument_group(description = "Parameters for Newton-Raphson iterations to estimate saturation time")
    nrp.add_argument("-A","--NR_alpha",type=float,default=1.)
    nrp.add_argument("-P","--NR_precision",type=float,default=1e-10)
    nrp.add_argument("-M","--NR_maxsteps",type=int,default=10000)
    return p


def PoissonSeedingVectors(m,n,cutoff = 1e-100,diff = False):
    px = np.zeros((len(n),len(m)))
    if diff:
        dpx = np.zeros((len(n),len(m)))
    for i in range(len(n)):
        if n[i] > 0:
            px[i] = poisson.pmf(m,n[i])
            px[i,px[i,:]<cutoff] = 0.
            px[i,-1] += (1. - np.sum(px[i]))
            if diff:
                dpx[i] = (m/n[i] - 1.)*px[i]
        else:
            px[i,0] = 1.
            if diff:
                dpx[i,1] = 1.
    if diff:
        return px,dpx
    else:
        return px


class MicrobialStrain(object):
    '''
    Stores all characteristics of a microbial strain
    mostly used to have always a consitent set of parameters
    (a strain has to bi initialized with each parameter)
    and also to check upon correct values when initializing
    or changing these parameters later
    
    so far, we implemented:
        * growth rate
        * yield factor
        * death rate
        
    '''
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
            super(MicrobialStrain,self).__setattr__(key,value)

            
class Environment(object):
    '''
    Class to store environmental parameters
    '''
    def __init__(self,substrate = 1e4,dilution = 1.,mixingtime = 10.,numdroplets = 1000):
        self.substrate   = substrate
        self.dilution    = dilution
        self.mixingtime  = mixingtime
        if not numdroplets is None:
            self.numdroplets = numdroplets
            self.__usenumdroplets = True
        else:
            self.__usenumdroplets = False
        
    def __getattr__(self,key):
        if key == "substrate":
            return self.__substrate
        elif key == "dilution":
            return self.__dilution
        elif key == "mixingtime":
            return self.__mixingtime
        elif key == "numdroplets":
            if self.__usenumdroplets:
                return self.__numdroplets
            else:
                return None
    
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
        elif key == "numdroplets" and self.__usenumdroplets:
            try:
                self.__numdroplets = int(value)
            except:
                raise ValueError
            if self.__numdroplets < 1:
                self.__numdroplets = 1
        else:
            super(Environment,self).__setattr__(key,value)
    
    def getParams():
        return {"substrate":self.substrate,"dilution":self.dilution,"mixingtime":self.mixingtime,"numdroplets":self.numdroplets}


class GrowthDynamics(object):
    def __init__(self,NR_alpha = 1.,NR_precision = 1e-14, NR_maxsteps = 10000,numstrains = None,**kwargs):
        
        if not numstrains is None:
            defaultlength = numstrains
        else:
            defaultlength = 1
        growthrates  = kwargs.get("growthrates",np.ones(defaultlength))
        yieldfactors = kwargs.get("yieldfactors",np.ones(defaultlength))
        assert len(growthrates) == len(yieldfactors)
        defaultlength = len(growthrates)
        if hasattr(kwargs,"deathrates"):
            deathrates = kwargs.get("deathrates")
            assert len(growthrates) == len(deathrates)
            self.__usedeathreates = True
        else:
            self.__usedeathreates = False
            deathrates = np.zeros(defaultlength)

        self.strains = list()
        for a,y,d in zip(growthrates,yieldfactors,deathrates):
            self.addStrain(growthrate = a,yieldfactor = y,deathrate = d)
        
        self.env = Environment( dilution = kwargs.get("dilution",1e-4),
                                mixingtime = kwargs.get("mixingtime",12),
                                substrate = kwargs.get("substrateconcentration",1e4),
                                numdroplets = kwargs.get("numdroplets") )
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
            # check if initialcells can be cast to array of floats
            ret_ic = np.array(initialcells,dtype=float)
        except:
            # fall back to (1,...,1) if cast does not work
            ret_ic = np.ones(self.numstrains,dype=float)
        if len(ret_ic) < self.numstrains:
            # fill up list of initial conditions with zeros
            ret_ic = np.concatenate((ret_ic,np.zeros(self.numstrains - len(ret_ic),dtype=float)))
        elif len(ret_ic) > self.numstrains:
            # or crop list if it is too long
            ret_ic = ret_ic[:self.numstrains]
        return ret_ic
        

    def Growth(self,initialcells = None):
        ic  = self.checkInitialCells(initialcells) # generate list with same dimensions as number of microbial strains
        ttd = self.__getTimeToDepletion(ic)          # time to depletion
        return self.env.dilution * ic * np.exp(self.growthrates * ttd - self.deathrates * self.env.mixingtime)


    def getGrowthVector(self,size):
        if isinstance(size,(int,np.int,np.int32,np.int64)):
            g = np.zeros(size)
            for i in np.arange(size):
                g[i] = self.Growth([i])
        elif isinstance(size,np.ndarray):
            g = np.zeros(len(size))
            i = 0
            for m in size:
                g[i] = self.Growth([m])
                i += 1
        return g

    def getGrowthMatrix(self,size):
        if isinstance(size,int):
            m = np.arange(size)
            g = np.zeros((size,size,2))
            for i in m:
                for j in m:
                    g[i,j] = self.Growth(initialcells = np.array([i,j]))
            return g[:,:,0],g[:,:,1]        
        elif isinstance(size,np.ndarray):
            if isinstance(size[0],np.ndarray):
                if len(size) >= 2:
                    m0 = size[0]
                    m1 = size[1]
                    g = np.zeros((len(m0),len(m1),2))
                    for i in range(len(m0)):
                        for j in range(len(m1)):
                            g[i,j] = self.Growth(initialcells = np.array([m0[i],m1[j]]))
                    return g[:,:,0],g[:,:,1]        
            else:
                m = size
                g = np.zeros(len(m))
                for i in range(len(m)):
                    g[i] = self.Growth(initialcells = np.array([m[i]]))[0]
                return g
        elif isinstance(size,(list,tuple)):
            if (len(size) >= 2) and isinstance(size[0],np.ndarray):
                m0 = size[0]
                m1 = size[1]
                g = np.zeros((len(m0),len(m1),2))
                for i in range(len(m0)):
                    for j in range(len(m1)):
                        g[i,j] = self.Growth(initialcells = np.array([m0[i],m1[j]]))
                return g[:,:,0],g[:,:,1]
                
    
    def getGrowthMultipleStrains(self,size,nstrains=2):
        g = [np.zeros(np.repeat(size,nstrains)) for i in range(nstrains)]
        for ic in itertools.product(range(size),repeat=nstrains):
            tmpgrowth = self.Growth(ic)
            for i in range(nstrains):
                g[i][ic] = tmpgrowth[i]
        return g
    
        
    def getSingleStrainFixedPoints(self):
        t = 1./self.growthrates * np.log(1./self.env.dilution)
        n = np.array([ self.yieldfactors[i] if t[i] <= self.env.mixingtime else 0. for i in range(self.numstrains)])
        if self.env.dilution < 1.:
            return self.env.dilution / (1. - self.env.dilution) * self.env.substrate * n
        else:
            return None

    def getSingleStrainFixedPointsApproximate(self):
        # approximate Poisson seeding with single strains.
        param = self.getSingleStrainFixedPoints()
        n = param - np.exp(-param+1)
        n[param<1] = 0
        return n
    
    def getSingleStrainFixedPointsPoissonSeeding(self,size=100):
        n = self.getSingleStrainFixedPointsApproximate()
        dn = 1.
        m = np.arange(size)
        for i in range(self.numstrains):
            if n[i] > 0.:
                growthi = np.zeros(size)
                for j in range(size):
                    ic = np.zeros(self.numstrains)
                    ic[i] = j
                    growthi[j] = self.Growth(ic)[i]
                step = 0
                dn = 1.
                # Newton-Raphson iteration to determine fixed point of non-linear equation (with Poisson seeding)
                while (dn/n[i])**2 > self.NR['precision2']:
                    px,dpx = PoissonSeedingVectors(m,np.array([n[i]]),diff=True)
                    dn = (np.dot(px[0],growthi) - n[i])/(np.dot(dpx[0],growthi) - 1.)
                    n[i] -= self.NR['alpha'] * dn
                    step += 1
                    if step > self.NR['maxsteps']:
                        break
                    if n[i] <= 0:
                        n[i] = 0
                        break
        return n
    
    def getTimeToDepletionMatrix(self,size):
        m = np.arange(size)
        t = np.zeros((size,size))
        for i in m:
            for j in m:
                t[i,j] = self.__getTimeToDepletion(initialcells = np.array([i,j]))
        return t

    def getApproximateGamma(self,initialcells):
        ic = self.checkInitialCells(initialcells)
        
        if (int(ic[0]) == 1) and (int(ic[1]) == 1):
            if self.growthrates[0] > self.growthrates[1]:
                invading = 1
            else:
                invading = 0
        elif int(ic[0]) == 1:
            invading = 0
        elif int(ic[1]) == 1:
            invading = 1
        else:
            raise ValueError
        
        noninvading = 1-invading
        gamma = np.zeros(2)
        a = self.growthrates[invading]/self.growthrates[noninvading]
        
        if a < 1:
            gamma[noninvading] = 1. - np.pow(self.env.substrate*self.yieldfactors[noninvading]/ic[noninvading],a)/(self.env.substrate*self.yieldfactor[invading])
        elif a==1:
            gamma[noninvading] = ic[noninvading]/(ic[noninvading] + 1)
        else:
            gamma[noninvading] = np.pow(ic[noninvading],a/(a-1.))/(self.env.substrate*self.yieldfactors[noninvading])*np.pow(self.yieldfactors[invading]/self.yieldfactors[noninvading],(a+1)/a)*np.pow(self.env.substrate*self.yieldfactors[noninvading]*np.pow(ic[noninvading],-a/(a-1))-1,1/a)
        
        gamma[invading] = 1 - gamma[noninvading]
        
        return gamma
            
            
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
            super(GrowthDynamics,self).__setattr__(key,value)
        
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
            super(StochasticGrowthDynamics,self).__getattr__(self,key)




class TimeIntegrator(object):
    # General forward integration of dynamics with Runge-Kutta method of 4th order
    # allows definition of multiple endconditions, currently implemented maximum time and one of the populations reaching zero
    def __init__(self,step = 1e-3,requiredpositive = True,initialconditions = None,dynamics = None,globaltime = 0,**kwargs):
        self.__step = step
        self.__requiredpositive = requiredpositive

        self.params = kwargs.get('params',None)
        
        if initialconditions is None:   raise ValueError
        else:                           self.x   = np.array(initialconditions)
        if dynamics is None:            raise NotImplementedError
        else:                           self.dyn = dynamics
        
        assert len(self.x) == len(self.dyn(0,self.x,self.params)), "Dimensions of initial conditions and dynamics do not match"
            
        self.__globaltime = globaltime        
        self.__EndConditions = list()
        self.__triggeredEndConditions = list()
        
    def RungeKutta4(self,xx,tt):
        # 4th order Runge-Kutta integration scheme
        k1 = self.__step * self.dyn( tt               , xx      , self.params )
        k2 = self.__step * self.dyn( tt+self.__step/2., xx+k1/2., self.params )
        k3 = self.__step * self.dyn( tt+self.__step/2., xx+k2/2., self.params )
        k4 = self.__step * self.dyn( tt+self.__step   , xx+k3   , self.params )
        return xx + (k1+2*k2+2*k3+k4)/6.

    def IntegrationStep(self,time):
        self.__triggeredEndConditions = list()
        t = 0
        while t <= time:
            self.x = self.RungeKutta4(self.x,self.__globaltime + t)
            if self.__requiredpositive:
                self.x[self.x<=0]=0
            t += self.__step
        self.__globaltime += time
    
    def IntegrateToZero(self,index):
        t = 0
        while self.x[index] > 0:
            self.x = self.RungeKutta4(self.x,self.__globaltime + t)
            if self.__requiredpositive:
                self.x[self.x<=0]=0
            t += self.__step
        self.__globaltime += t
        self.__triggeredEndConditions = list(["reachzero",index])
        
    def ResetInitialConditions(self,initialconditions,globaltime = 0):
        self.x = np.array(initialconditions,dtype=np.float64)
        assert len(self.x) == len(self.dyn(0,self.x,self.params)), "Dimensions of initial conditions and dynamics do not match"
        self.__globaltime = globaltime
        self.__triggeredEndConditions = list()

    def SetEndConditionMaxTime(self,maxtime):
        if float(maxtime) >= 0:
            self.__EndConditions.append(["maxtime",float(maxtime)])
        else:
            raise ValueError
    
    def SetEndConditionReachZero(self,populationindex):
        if len(self.x) <= populationindex:
            raise IndexError
        self.__EndConditions.append(["reachzero",populationindex])
        
    def SetEndCondition(self,condition,value):
        if str(condition).lower() == "maxtime":
            if float(value) >= 0:
                self.__EndConditions.append(["maxtime",float(value)])
        elif str(condition).lower() == "reachzero":
            if len(self.x) > int(value):
                self.__EndConditions.append(["reachzero",int(value)])
        else:
            raise NotImplementedError

    def HasEnded(self):
        terminateInteration = False
        if np.isnan(self.x).any():
            terminateInteration = True
        else:
            for ec in self.__EndConditions:
                if ec[0] == "maxtime":
                    if ec[1] < self.__globaltime:
                        terminateInteration = True
                        self.__triggeredEndConditions.append(ec)
                elif ec[0] == "reachzero":
                    if self.x[ec[1]] <= 0.:
                        terminateInteration = True
                        self.__triggeredEndConditions.append(ec)
                else:
                    raise NotImplementedError
        return terminateInteration
    
    def IntegrateToEndConditions(self):
        if self.CountEndConditions > 0:
            while not self.HasEnded():
                self.x = self.RungeKutta4(self.x,self.__globaltime)
                if self.__requiredpositive:
                    self.x[self.x<=0]=0
                self.__globaltime += self.__step
        else:
            raise NotImplementedError

    def __str__(self):
        return (" ".join(["{:14.6e}"]*len(self.x))).format(*self.x)
    
    
    def __getattr__(self,key):
        if key == "CountEndConditions":
            return len(self.__EndConditions)
        elif key == "populations":
            return self.x
        elif key == "time":
            return self.__globaltime
        else:
            super(TimeIntegrator,self).__getattr__(self,key)
    
    def __getitem__(self,key):
        if int(key) < len(self.x):
            return self.x[int(key)]


class GrowthDynamicsPublicGoods(GrowthDynamics):
    def __init__(self,numstrains = None,**kwargs):
        
        if kwargs.get("mixingtime") is None:
            kwargs["mixingtime"] = 12.
        
        super(GrowthDynamicsPublicGoods,self).__init__(self,numstrains = numstrains,**kwargs)
        
        if kwargs.get("polynomialinteraction",True):
            # polynomial coefficients for interaction with public good
            # a = Sum_n a_n G^n
            # y = Sum_n y_n G^n
            self.__PGInteractionGrowthRates = np.array(kwargs.get("pginteractiongrowthrates",np.zeros(self.numstrains)),dtype=np.float64)
            self.__PGInteractionYieldFactor = np.array(kwargs.get("pginteractionyieldfactor",np.zeros(self.numstrains)),dtype=np.float64)
            self.__PGGrowthRatesOrder = len(self.__PGInteractionGrowthRates)/self.numstrains
            self.__PGYieldFactorOrder = len(self.__PGInteractionYieldFactor)/self.numstrains
            self.__PGInteractionGrowthRates = np.reshape(self.__PGInteractionGrowthRates,(self.numstrains,self.__PGGrowthRatesOrder))
            self.__PGInteractionYieldFactor = np.reshape(self.__PGInteractionYieldFactor,(self.numstrains,self.__PGYieldFactorOrder))
            
            self.GR = self.PolynomialGrowthRates
            self.YF = self.PolynomialYieldFactors
            
        else:
            # exponential interaction with public good
            # a = a * (A + (1-2A) EXP(-eps G))
            # y = y * (B + (1-2B) EXP(-delta G))
            # coefficients eps, delta obtained from commandline parameters
            # if commandline parameters have 2n numbers, second set determines A, B 
            #         exponentially decreasing for 0, exponentially increasing for 1 (default: vectors of 0)
            self.__PGInteractionGrowthRates = np.array(kwargs.get("pginteractiongrowthrates",np.zeros(self.numstrains)),dtype=np.float64)
            self.__PGInteractionYieldFactor = np.array(kwargs.get("pginteractionyieldfactor",np.zeros(self.numstrains)),dtype=np.float64)
            
            if len(self.__PGInteractionGrowthRates) == self.numstrains:
                self.__PGInteractionGrowthRates = np.array([self.__PGInteractionGrowthRates,np.zeros(self.numstrains)])
            elif len(self.__PGInteractionGrowthRates) == 2 * self.numstrains:
                self.__PGInteractionGrowthRates = np.reshape(self.__PGInteractionGrowthRates,(self.numstrains,2))
            else:
                raise ValueError
            
            if len(self.__PGInteractionYieldFactor) == self.numstrains:
                self.__PGInteractionYieldFactor = np.array([self.__PGInteractionYieldFactor,np.zeros(self.numstrains)])
            elif len(self.__PGInteractionYieldFactor) == 2 * self.numstrains:
                self.__PGInteractionYieldFactor = np.reshape(self.__PGInteractionYieldFactor,(self.numstrains,2))
            else:
                raise ValueError
            
            self.GR = self.ExponentialGrowthRates
            self.YF = self.ExponentialYieldFactors
            
            
        
        self.__PGProduction = np.array(kwargs.get("pgproduction",np.zeros(self.numstrains)),dtype=np.float64)
        assert len(self.__PGProduction) == self.numstrains, "production of public goods does not match number of strains"
        assert sum(self.__PGProduction) > 0, "no public goods produced"

        self.dyn = TimeIntegrator(dynamics = self.PGdyn,initialconditions = np.ones(self.numstrains+2),params = None)
        self.dyn.SetEndCondition("maxtime",self.env.mixingtime)
        self.dyn.SetEndCondition("reachzero",self.numstrains)
        
        self.__onlypositivecoefficients = kwargs.get("onlypositivecoefficients",True)
        
    
    # dynamics for all strains, then substrate, then public good
    def PGdyn(self,t,x,params):
        # public good can influence growth rates and yield
        a = self.GR(x)
        y = self.YF(x)
        return np.concatenate([a*x[:-2],np.array([-np.sum(a/y*x[:-2]),np.sum(self.__PGProduction*x[:-2])])]) # cellcounts, substrate, pg
    
    
    def Growth(self,initialcells = None):
        ic = self.checkInitialCells(initialcells)
        ic = np.concatenate([ic,np.array([self.env.substrate,0])]) # cellcounts, substrate, pg
        self.dyn.ResetInitialConditions(ic)
        self.dyn.IntegrateToEndConditions()
        return self.dyn.populations[:-2] # return only cell count
    
    # polynomial dependence on public good concentration
    def PolynomialGrowthRates(self,populations):
        da = np.zeros(self.numstrains)
        # start with highest order and go back to first
        for i in range(1,self.__PGGrowthRatesOrder+1):
            da += self.__PGInteractionGrowthRates[:,-i]
            da *= populations[-1]
        a = self.growthrates + da
        if self.__onlypositivecoefficients:
            a[a<0] = 0
        return a
    
    # polynomial dependence on public good concentration
    def PolynomialYieldFactors(self,populations):
        dy = np.zeros(self.numstrains)
        # start with highest order and go back to first
        for i in range(1,self.__PGYieldFactorOrder+1):
            dy += self.__PGInteractionYieldFactor[:,-i]
            dy *= populations[-1]
        y = self.yieldfactors + dy
        if self.__onlypositivecoefficients:
            y[y<1e-300] = 1e-300
        return y
    
    def ExponentialGrowthRates(self,populations):
        return self.growthrates * (self.__PGInteractionGrowthRates[1,:] + (1-2*self.__PGInteractionGrowthRates[1,:])*np.exp(self.__PGInteractionGrowthRates[0,:] * populations[-1]))
    
    def ExponentialYieldFactors(self,populations):
        return self.yieldfactors * (self.__PGInteractionYieldFactor[1,:] + (1-2*self.__PGInteractionYieldFactor[1,:])*np.exp(self.__PGInteractionYieldFactor[0,:] * populations[-1]))
    
                            
    





