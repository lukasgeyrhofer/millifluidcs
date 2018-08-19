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


  Lukas Geyrhofer, l.geyrhofer@technion.ac.il, 2016-2018

'''

import numpy as np
import argparse
from scipy.stats import poisson
import itertools
import scipy.integrate as spint
import inspect
import sys

import pickle

def RungeKutta4(func,xx,tt,step):
  # 4th order Runge-Kutta integration scheme
  k1 = step * func( tt        , xx )
  k2 = step * func( tt+step/2., xx+k1/2. )
  k3 = step * func( tt+step/2., xx+k2/2. )
  k4 = step * func( tt+step   , xx+k3 )
  return xx + (k1+2*k2+2*k3+k4)/6.


def AddGrowthParameters(p,allparams = False,deathrates = False,numdroplets = False,dilution = False,
                        defaultgrowthrates = [2.,1.],defaultyieldfactors = [1.,2.],defaultdeathrates = [0.,0.],
                        defaultsubstrate = 1e4, defaultmixingtime = 24,defaultdilution = 2e-4, defaultnumdroplets = 1000):
    # Helper routine to generate all cmdline parameters for microbial growth
    gp = p.add_argument_group(description = "==== Parameters for growth in droplets ====")
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
    nrp = p.add_argument_group(description = "==== Parameters for Newton-Raphson iterations to estimate saturation time ====")
    nrp.add_argument("-A","--NR_alpha",type=float,default=1.)
    nrp.add_argument("-P","--NR_precision",type=float,default=1e-10)
    nrp.add_argument("-M","--NR_maxsteps",type=int,default=10000)
    return p


def PoissonSeedingVectors(m,n,cutoff = 1e-100,diff = False):
    if isinstance(n,(float,np.float,np.float64,int)):
        n = np.array([n],dtype=float)
    px = np.zeros((len(n),len(m)))
    if diff:
        dpx = np.zeros((len(n),len(m)))
    for i in range(len(n)):
        if n[i] > 0:
            px[i] = poisson.pmf(m,n[i])
            px[i,px[i,:]<cutoff] = 0.
            px[i] /= np.sum(px[i]) # normalize
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


def PointSeedingVectors(m,n):
    if isinstance(n,(float,np.float,np.float64,int)):
        n = np.array([n],dtype=float)
    px = np.zeros((len(n),len(m)))
    for i in range(len(n)):
        if n[i] > 0:
            idx0 = ((m - n[i])**2).argmin()
            try:
                if m[idx0] < n[i] and m[idx0+1] > n[i]:
                    p[idx0] = (m[idx0+1] - n[i])/(m[idx0+1] - m[idx0])
                    p[idx0+1] = 1 - p[idx0]
                if m[idx0-1] < n[i] and m[idx0] > n[i]:
                    p[idx0-1] = (m[idx0] - n[i])/(m[idx0] - m[idx0-1])
                    p[idx0] = 1 - p[idx0-1]
            except:
                pass
    return px


def getInoculumAxes(**kwargs):
    newcoord  = kwargs.get("newcoordinates",False)
    nmax      = kwargs.get("maxInoculum",50)
    nstep     = kwargs.get("stepInoculum",2)
    nlist     = np.arange(start = 0,stop = nmax + .5*nstep,step = nstep)
    if newcoord:
        xstep = kwargs.get("stepFraction",0.05)
        xlist = np.arange(start = 0,stop = 1 + .5*xstep,step = xstep)
        return nlist,xlist
    else:
        return nlist,nlist

def getAbsoluteInoculumNumbers(coordinate,newcoordinates = False):
    if newcoordinates:
        return coordinate[0] * coordinate[1], coordinate[0] * (1 - coordinate[1])
    else:
        return coordinate[0],coordinate[1]

def SeedingAverage(matrix,coordinates,axis1 = None, axis2 = None, mask = None, replaceNAN = True):
    dim = matrix.shape
    if axis1 is None:   axis1 = np.arange(dim[0])
    if axis2 is None:   axis2 = np.arange(dim[1])
        
    p1 = PoissonSeedingVectors(axis1,[coordinates[0]])[0]
    p2 = PoissonSeedingVectors(axis2,[coordinates[1]])[0]
    
    matrix0 = matrix[:,:]
    if not mask is None:    matrix0[np.logical_not(mask)] = 0
    if replaceNAN:          matrix0 = np.nan_to_num(matrix0)
    
    return np.dot(p2,np.dot(p1,matrix0))
    



def AssignGrowthDynamics(**kwargs):
    # pick GrowthDynamics class from below via argument string
    # convert all values of the dict-entry 'ParameterList' into entries of kwargs itself

    def AddEntry(d,key,val):
        tmp = dict()
        if not key is None:
            if len(val) == 1:
                tmp[key] = val[0]
            elif len(val) > 1:
                tmp[key] = np.array(val)
        tmp.update(d)
        return tmp

    def MakeDictFromParameterList(params):
        p = dict()
        curkey = None
        curvalue = list()
        for entry in params:
            try:
                v = float(entry)
                curvalue.append(v)
            except:
                p = AddEntry(p,curkey,curvalue)
                curvalue = list()
                curkey = entry

        p = AddEntry(p,curkey,curvalue)
        return p
    
    GrowthDynamics = kwargs.get('GrowthDynamics','')
    
    # generate dict from parameters and update with kwargs
    params = MakeDictFromParameterList(kwargs.get('ParameterList',[]))
    params.update(kwargs)
    
    # get rid of original description for these parameters,
    # such that they are not passed twice in different form
    if 'ParameterList' in params.keys():
        del params['ParameterList']
    
    for name,dyn in inspect.getmembers(sys.modules['growthclasses'],inspect.isclass):
        if name == 'GrowthDynamics' + GrowthDynamics.strip():
            return dyn(**params)
    
    # did not find GrowthDynamics
    raise NotImplementedError("'GrowthDynamics{}' not yet implemented.".format(GrowthDynamics.strip()))


def AddGrowthDynamicsArguments(p):
    c = [name[14:] for name,obj in inspect.getmembers(sys.modules['growthclasses'],inspect.isclass) if name[:14] == 'GrowthDynamics']
    pgd = p.add_argument_group(description = "==== GrowthDynamics ====")
    pgd.add_argument("-G","--GrowthDynamics",choices = c, default = '')
    pgd.add_argument("-P","--ParameterList",nargs="*",default = [])
    return p


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
        self.__growthrate  = growthrate
        self.__yieldfactor = yieldfactor
        self.__deathrate   = deathrate
        
        self.__growing   = True
    
    def StopGrowth(self):
        self.__growing = False
    def AllowGrowth(self):
        self.__growing = True
    
    def __getattr__(self,key):
        if key == "growthrate":
            if self.__growing:
                return self.__growthrate
            else:
                return 0.
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
    def __init__(self,numstrains = None,**kwargs):
        
        if not numstrains is None:
            defaultlength = numstrains
        else:
            defaultlength = 1
        growthrates  = kwargs.get("growthrates",np.ones(defaultlength))
        yieldfactors = kwargs.get("yieldfactors",np.ones(defaultlength))
        assert len(growthrates) == len(yieldfactors),"All strains need growthrate and yield defined"
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
        
        self.env = Environment( dilution    = kwargs.get("dilution",               1.),
                                mixingtime  = kwargs.get("mixingtime",             24),
                                substrate   = kwargs.get("substrateconcentration", 1e4),
                                numdroplets = kwargs.get("numdroplets") )
        
        self.NR  = {'alpha' :       float(kwargs.get("NR_alpha",     1. )),
                    'precision2' :  float(kwargs.get("NR_precision", 1e-14 ))**2,
                    'maxsteps' :    int(  kwargs.get("NR_maxsteps",  1000 ))}

        
        self.__growthmatrix      = None
        self.__growthmatrixgridX = None
        self.__growthmatrixgridY = None
        
        self.__kwargs_for_pickle = kwargs
        
        self.__hasximatrix = False
        self.__restrictFractionalPopulation = kwargs.get("RestrictFractionalPopulation",True)
        
        self.__trajectorytimestep = kwargs.get('TimeIntegratorStep',1e-3) * kwargs.get('TimeIntegratorOutput',10)
        self.__params = dict()
        
    
    def addStrain(self,growthrate = 1.,yieldfactor = 1.,deathrate = 0):
        self.strains.append(MicrobialStrain(growthrate = growthrate, yieldfactor = yieldfactor, deathrate = deathrate))
    
    def delLastStrain(self):
        return self.strains.pop()
    
    def AllowGrowth(self):
        for i in range(len(self.strains)):
            self.strains[i].AllowGrowth()
    def StopGrowth(self,strain = None):
        if strain is None:
            for i in range(len(self.strains)):
                self.strains[i].StopGrowth()
        else:
            if strain < len(self.strains):
                self.strains[strain].StopGrowth()
    
    def getTimeToDepletion(self,initialcells):
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
        if initialcells is None:
            ret_ic = np.ones(self.numstrains,dtype = float)
        else:
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


    def Trajectory(self,initialconditions,TimeOutput = False):
        ic = self.checkInitialCells(initialconditions)
        tdepl = self.getTimeToDepletion(ic)
        t = np.arange(start = self.__trajectorytimestep, stop = self.env.mixingtime, step = self.__trajectorytimestep)
        x = list([ic])
        for tcur in t:
            if tcur < tdepl:
                x.append(ic * np.exp(self.growthrates * tcur))
            else:
                x.append(ic * np.exp(self.growthrates * tdepl))
        xx = np.vstack(x)
        if TimeOutput:
            tt = np.array([np.concatenate([[0],t])]).T
            return np.concatenate([tt,xx],axis = 1)
        else:
            return xx
        

    def Growth(self,initialcells = None):
        ic  = self.checkInitialCells(initialcells) # generate list with same dimensions as number of microbial strains
        ttd = self.getTimeToDepletion(ic)          # time to depletion
        g   = self.env.dilution * ic * np.exp(self.growthrates * ttd - self.deathrates * self.env.mixingtime)
        return g


    def getGrowthVector(self,size,strainID = 0):
        if isinstance(size,(int,np.int,np.int32,np.int64)):
            g = np.zeros(size)
            m = np.zeros(self.numstrains)
            for i in np.arange(size):
                m[strainID] = i
                g[i] = self.Growth(m)[strainID]
        elif isinstance(size,np.ndarray):
            g = np.zeros(len(size))
            m = np.zeros(self.numstrains)
            i = 0
            for j in size:
                m[strainID] = j
                g[i] = self.Growth(m)[strainID]
                i += 1
        return g


    def ComputeGrowthMatrix(self,size,step=1):
        if isinstance(size,int):
            self.__growthmatrixgridX = np.arange(start = 0, stop = size, step = step)
            self.__growthmatrixgridY = np.arange(start = 0, stop = size, step = step)
        elif isinstance(size,(list,tuple,np.ndarray)):
            if isinstance(size[0],int):
                self.__growthmatrixgridX = np.arange(start = 0,stop = size[0],step = step)
            elif isinstance(size[0],(list,tuple,np.ndarray)):
                self.__growthmatrixgridX = size[0]
            else:
                raise ValueError("size argument can only be int or (list/tuple of int)")

            if len(size) >= 2:
                if isinstance(size[1],int):
                    self.__growthmatrixgridY = np.arange(start = 0,steop = size[1],step = step)
                elif isinstance(size[1],(list,tuple,np.ndarray)):
                    self.__growthmatrixgridY = size[1]
                else:
                    raise ValueError("size argument can only be int or (list/tuple of int)")
            else:
                self.__growthmatrixgridY = self.__growthmatrixgridX[:]
                
        else:
            raise ValueError("size argument does not fit")

        self.__growthmatrix = np.zeros((len(self.__growthmatrixgridX),len(self.__growthmatrixgridY),2))
        for i,n1 in enumerate(self.__growthmatrixgridX):
            for j,n2 in enumerate(self.__growthmatrixgridY):
                self.__growthmatrix[i,j] = self.Growth(initialcells = np.array([n1,n2]))

    
    def getGrowthMatrix(self,size,step=1):
        # backwards compatibility
        self.ComputeGrowthMatrix(size,step)
        return self.__growthmatrix
    
    def hasGrowthMatrix(self):
        return not (self.__growthmatrix is None)
    
    def ExtendGrowthMatrix(self,size,step=1):
        if isinstance(size,int):
            if size > self.__growthmatrixgridX[-1]:
                new_growthmatrixgridX = np.concatenate((self.__growthmatrixgridX,np.arange(start = self.__growthmatrixgridX[-1]+step,stop = size,step = step)))
            else:
                new_growthmatrixgridX = self.__growthmatrixgridX
            if size > self.__growthmatrixgridY[-1]:
                new_growthmatrixgridY = np.concatenate((self.__growthmatrixgridY,np.arange(start = self.__growthmatrixgridY[-1]+step,stop = size,step = step)))
            else:
                new_growthmatrixgridY = self.__growthmatrixgridY
        else:
            raise NotImplementedError

        g = np.zeros((len(new_growthmatrixgridX),len(new_growthmatrixgridY),2))
        for i in range(len(new_growthmatrixgridX)):
            x = new_growthmatrixgridX[i]
            for j in range(len(new_growthmatrixgridY)):
                y = new_growthmatrixgridY[j]
                if (x in self.__growthmatrixgridX) and (y in self.__growthmatrixgridX):
                    g[i,j] = self.__growthmatrix[i,j]
                else:
                    g[i,j] = self.Growth(initialcells = np.array([x,y]))
                    
        self.__growthmatrixgridX = new_growthmatrixgridX[:]
        self.__growthmatrixgridY = new_growthmatrixgridY[:]
        self.__growthmatrix      = g[:,:,:]
    
    
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
                t[i,j] = self.getTimeToDepletion(initialcells = np.array([i,j]))
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


    def setGrowthMatrixValues(self,threshold,newvalue = 0, operation = 'below'):
        if self.hasGrowthMatrix():
            if operation.lower() == 'below':
                self.__growthmatrix[self.__growthmatrix <  threshold] = newvalue
            elif operation.lower() == 'above':
                self.__growthmatrix[self.__growthmatrix >  threshold] = newvalue
            elif operation.lower() == 'equal':
                self.__growthmatrix[self.__growthmatrix == threshold] = newvalue
            else:
                raise NotImplementedError


    def GetXi(self,initialconditions):
        if self.hasGrowthMatrix():
            
            # reverse compute expansion factor xi from total population size N,
            # using the analytic solution N(t) = n xi SUM_j[x_j xi(t)^da_j]
            
            ic   = self.checkInitialCells(initialconditions)
            idx0 = ((self.growthmatrixgrid[0] - ic[0])**2).argmin()
            idx1 = ((self.growthmatrixgrid[1] - ic[1])**2).argmin()
            n    = np.sum(ic)
            
            if n>0:
                x       = ic/n
                nfin    = np.sum(self.growthmatrix[idx0,idx1,:])
                xi0     = nfin/n
                da      = self.growthrates/np.mean(self.growthrates) - 1.
                xi      = xi0
                xi_last = 0
                i       = 0
                while ((xi_last-xi)**2) > self.NR['precision2'] * (xi**2):
                    # store value to measure convergence
                    xi_last = xi
                    
                    # Newton-Raphson iteration to refine solution
                    xida   = np.power(xi, da)
                    Sxi    = xi * np.dot(x, xida)
                    Sxddxi = np.dot(x * (1. + da), xida)
                    xi    -= self.NR['alpha']*(Sxi - xi0)/Sxddxi
                    
                    # should not iterate infinitely
                    i     += 1
                    if i > self.NR['maxsteps']:
                        raise ValueError
                    
                return xi
            else:
                return 0.
        else:
            raise ValueError

    def ComputeXiMatrix(self):
        if self.hasGrowthMatrix:
            self.__xi = np.empty((len(self.growthmatrixgrid[0]),len(self.growthmatrixgrid[1])))
            for i,n1 in enumerate(self.growthmatrixgrid[0]):
                for j,n2 in enumerate(self.growthmatrixgrid[1]):
                    self.__xi[i,j] = self.GetXi([n1,n2])
            self.__hasximatrix = True
        
    def GetXiMatrix(self):
        if not self.__hasximatrix:
            self.ComputeXiMatrix()
        return self.__xi

    def __getattr__(self,key):
        if key == "numstrains":
            return len(self.strains)
        elif key == "growthrates":
            return np.array([self.strains[i].growthrate for i in range(self.numstrains)])
        elif key == "yieldfactors":
            return np.array([self.strains[i].yieldfactor for i in range(self.numstrains)])
        elif key == "deathrates":
            return np.array([self.strains[i].deathrate for i in range(self.numstrains)])
        elif key == "growthmatrix":
            if self.__growthmatrix is None:
                raise ValueError("Growthmatrix not yet computed")
            else:
                tmp = self.__growthmatrix
                if self.__restrictFractionalPopulation:
                    tmp[tmp<1] = 0
                return tmp
        elif key == "growthmatrixgrid":
            if (self.__growthmatrixgridX is None) or (self.__growthmatrixgridY is None):
                raise ValueError("Growthmatrix not yet computed")
            else:
                return (self.__growthmatrixgridX,self.__growthmatrixgridY)
        #else:
            #super(GrowthDynamics,self).__getattr__(key,value)


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
            
    
    def arraystring(self,x):
        return "[" + ", ".join(["{:.4f}".format(a) for a in x]) + "]"

    def ParameterString(self):
        r  = '\n'
        s  = "*** microbial strains ***" +r
        s += "  growthrates " + self.arraystring(self.growthrates) +r
        s += "  yield       " + self.arraystring(self.yieldfactors) +r+r
        s += "*** environment ***" +r
        s += "  mixingtime  " + str(self.env.mixingtime) +r
        s += "  substrate   " + str(self.env.substrate) +r
        if self.env.dilution < 1:
            s += "  dilution    " + str(self.env.dilution) +r
        return s

    def __str__(self):
        return self.ParameterString()

    # pickle functions for saving and loading object from file
    def __getstate__(self):
        return [self.__kwargs_for_pickle,self.__growthmatrix,(self.__growthmatrixgridX,self.__growthmatrixgridY)]
    
    def __setstate__(self,state):
        # python 3 keys in dictionaries are <type 'bytes'>, but in order to use **kwargs, we need <type 'str'>
        # need to check conversion
        kwargs = dict()
        for k,v in state[0].items():
            if isinstance(k,str):   kwargs[k] = v
            else:                   kwargs[k.decode('utf-8')] = v
        self.__init__(**kwargs)
        self.__growthmatrix = state[1]
        if isinstance(state[2],int):
            # backward compatibility
            self.__growthmatrixgridX = np.arange(state[2])
            self.__growthmatrixgridY = np.arange(state[2])
        else:
            # current implementation
            self.__growthmatrixgridX,self.__growthmatrixgridY = state[2]

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
    def __init__(self,step = 1e-3,requiredpositive = True,initialconditions = None,dynamics = None,globaltime = 0,outputstep = 100,**kwargs):
        self.__step = step
        self.__outputstep = int(outputstep)

        self.params = kwargs.get('params',None)
        
        self.have_start_values = False
        if not initialconditions is None:
            self.x   = np.array(initialconditions)
            self.have_start_values = True

        if dynamics is None:
            raise NotImplementedError
        else:
            self.dyn = dynamics
        
        #assert len(self.x) == len(self.dyn(0,self.x)), "Dimensions of initial conditions and dynamics do not match"
            
        self.__globaltime = globaltime        
        self.__EndConditions = list()
        self.__triggeredEndConditions = list()
        
        self.__extinctionthresholds = dict()
        self.__requiredpositive = requiredpositive
        self.__minimalpositivevalue = kwargs.get('MinimalPositiveValue',0)  # can introduce hard-cutoff, if values get too small
                                                                            # if requiredpositive is true, set everything below this threshold to 0
        
        self.__trajectory = None
    
    def RungeKutta4(self,xx,tt):
        # 4th order Runge-Kutta integration scheme
        k1 = self.__step * self.dyn( tt               , xx      )
        k2 = self.__step * self.dyn( tt+self.__step/2., xx+k1/2.)
        k3 = self.__step * self.dyn( tt+self.__step/2., xx+k2/2.)
        k4 = self.__step * self.dyn( tt+self.__step   , xx+k3   )
        ret = xx + (k1+2*k2+2*k3+k4)/6.
        if self.__requiredpositive:
            ret[ret < self.__minimalpositivevalue] = 0
        return ret
    
    
    def checkExtinction(self):
        if len(self.__extinctionthresholds) >= 1:
            for i in self.__extinctionthresholds.keys():
                if i < len(self.x):
                    if self.x[i] < self.__extinctionthresholds[i]:
                        self.x[i] = 0

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


    def IntegrationStep(self,time):
        if not self.have_start_values:
            raise ValueError
        
        self.__triggeredEndConditions = list()
        t = 0
        while t <= time:
            self.x = self.RungeKutta4(self.x,self.__globaltime + t)
            self.checkExtinction()
            t += self.__step
        self.__globaltime += t
        return self.__globaltime
    
    
    def IntegrateToZero(self,index):
        if not self.have_start_values:
            raise ValueError
        t = 0
        while self.x[index] > 0:
            self.x = self.RungeKutta4(self.x,self.__globaltime + t)
            self.checkExtinction()
            t += self.__step
        self.__globaltime += t
        self.__triggeredEndConditions = list(["reachzero",index])

    
    def IntegrateToEndConditions(self,store_trajectory = False):
        if not self.have_start_values:
            raise ValueError
        if store_trajectory:
            self.__trajectory = list()
        o = 0
        if self.CountEndConditions > 0:
            while not self.HasEnded():
                self.x = self.RungeKutta4(self.x,self.__globaltime)
                self.checkExtinction()
                self.__globaltime += self.__step
                if store_trajectory:
                    if o%self.__outputstep == 0:
                        self.__trajectory.append([self.__globaltime,self.x])
                o += 1
            return self.x
        else:
            raise NotImplementedError
        
    def ResetInitialConditions(self,initialconditions,globaltime = 0):
        self.x = np.array(initialconditions,dtype=np.float64)
        assert len(self.x) == len(self.dyn(0,self.x)), "Dimensions of initial conditions and dynamics do not match"
        self.__globaltime = globaltime
        self.__triggeredEndConditions = list()
        self.have_start_values = True

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

    
    def SetPopulation(self,index,value):
        if int(index) < len(self.x):
            self.x[int(index)] = float(value)
    
    def SetPopulationExtinctionThreshold(self,index,value):
        self.__extinctionthresholds[int(index)] = float(value)
    
    def __str__(self):
        return (" ".join(["{:14.6e}"]*len(self.x))).format(*self.x)
    
    def GetTrajectory(self,TimeOutput = False):
        if not self.__trajectory is None:
            if TimeOutput:
                return np.array([np.concatenate([np.array([t]),x]) for t,x in self.__trajectory])
            else:
                return np.array([x for t,x in self.__trajectory])
        else:
            raise ValueError
    
    def __getattr__(self,key):
        if key == "CountEndConditions":
            return len(self.__EndConditions)
        elif key == "populations":
            return self.x
        elif key == "time":
            return self.__globaltime
        #else:
            #super(TimeIntegrator,self).__getattr__(self,key)
    
    def __getitem__(self,key):
        if int(key) < len(self.x):
            return self.x[int(key)]

class GrowthDynamicsODE(GrowthDynamics):
    def __init__(self,numstrains = None, **kwargs):
        super(GrowthDynamicsODE,self).__init__(numstrains = numstrains,**kwargs)

        # mixingtime is crucial for integration
        # needs to be set by hand, if not provided as cmdline argument
        if self.env.mixingtime is None:
            self.env.mixingtime = 24.
        
        self.IntegrationMethod       = kwargs.get("IntegrationMethod",'OWNRK4')
        self.otherinitialconditions  = np.array([self.env.substrate],dtype=np.float64)
        self.TimeIntegratorOutput    = kwargs.get("TimeIntegratorOutput",10)
        self.TimeIntegratorStep      = kwargs.get("TimeIntegratorStep",1e-3)
        self.EmptySubstrateThreshold = 1e-3 * np.mean(self.yieldfactors)
        
        if self.IntegrationMethod.upper() == 'OWNRK4':
            # use TimeIntegrator class defined above
            self.integrator = TimeIntegrator(dynamics = self.dynamics,requiredpositive = True,**kwargs)
            self.Growth     = self.GrowthOwnRK4Integrator
            self.Trajectory = self.TrajectoryOwnRK4Integrator
            
            self.integrator.SetEndCondition("maxtime",self.env.mixingtime)
            for i in range(self.numstrains):
                self.integrator.SetPopulationExtinctionThreshold(i,1)

        elif self.IntegrationMethod.upper() == 'SCIPY':
            # initialize integration from 'Scipy.integrate' = 'spint'
            self.integrator = spint.ode(self.dynamics)
            self.integrator.set_integrator('vode', method = 'bdf', min_step = 1e-4, max_step = 1e-2)

            self.Trajectory = self.TrajectorySciPyIntegrator
            self.Growth     = self.GrowthSciPyIntegrator
        
        else:
            raise NotImplementedError
        
        self.AllowGrowth()
        
        
    # this function needs to be overwritten in all child-objects
    # see lines above, where this is set as the system of differential equations
    # child-objects only need to define this 'dynamics' to work
    def dynamics(self,t,x):
        a = self.growthrates
        if x[-1] <= 0:
            a = np.zeros(self.numstrains)
        return np.concatenate([
                    self.growthrates * x[:self.numstrains],
                    np.array([np.sum(-a * x[:self.numstrains]/self.yieldfactors)])
                ])

    # base growth function to use for time integrator dynamics
    def GrowthOwnRK4Integrator(self,initialcells = None):
        #ic = self.checkInitialCells(initialcells)
        #ic = np.concatenate([ic,self.otherinitialconditions])
        #self.integrator.ResetInitialConditions(ic)
        #final = self.integrator.IntegrateToEndConditions()
        #return final[:self.numstrains]
        
        # compute whole trajectory, only output final cell numbers
        tmp = self.TrajectoryOwnRK4Integrator(initialcells)
        return tmp[-1,:self.numstrains]
        


    # should also work for more complicated dynamics implemented in classes inherited from this one
    def TrajectoryOwnRK4Integrator(self,initialconditions,TimeOutput = False):
        # helper routine
        #def AddTimeToOutputVector(x,t,TimeOutput):
            #if TimeOutput:
                #return np.concatenate([np.array([t]),x])
            #else:
                #return x
        
        ## set initial conditions
        #initialconditions[:self.numstrains] = self.checkInitialCells(initialconditions[:self.numstrains])
        #if len(initialconditions) >= self.numstrains:
            #initialconditions = np.concatenate([initialconditions,self.otherinitialconditions])
        #self.integrator.ResetInitialConditions(initialconditions)
        
        ## generate first entry in output data
        #r = list()
        #r.append(AddTimeToOutputVector(self.integrator.x,0,TimeOutput))
        #i = 0
        #while not self.integrator.HasEnded():
            #t = self.integrator.IntegrationStep(self.TimeIntegratorStep)
            #i+=1
            #if i%self.TimeIntegratorOutput == 0:
                #r.append(AddTimeToOutputVector(self.integrator.x,t,TimeOutput))
        ## return output
        #return np.vstack(r)

        initialconditions[:self.numstrains] = self.checkInitialCells(initialconditions[:self.numstrains])
        if len(initialconditions) >= self.numstrains:
            initialconditions = np.concatenate([initialconditions,self.otherinitialconditions])
        self.integrator.ResetInitialConditions(initialconditions)
        self.integrator.IntegrateToEndConditions(store_trajectory = True)
        return self.integrator.GetTrajectory(TimeOutput)
        


    # compute a full trajectory, output the matrix of solutions
    def TrajectorySciPyIntegrator(self,initialcells,TimeOutput = False):

        # store output
        localtraj = []
        def solout(t,y):
            if TimeOutput:  localtraj.append(np.concatenate([[t],y]))
            else:           localtraj.append(y)
        
        # check number of cells
        # and append all other initialconditions and set initial conditions
        ic = self.checkInitialCells(initialcells[:self.numstrains])
        ic = np.concatenate([ic,self.otherinitialconditions])
        
        self.integrator.set_initial_value(ic,0)
        
        self.AllowGrowth()
        # integrate ODE
        while (self.integrator.t < self.env.mixingtime) and self.integrator.successful():
            self.integrator.integrate(self.integrator.t + self.TimeIntegratorStep)
            for strain in np.where(self.integrator.y[:self.numstrains] < 1)[0]:
                self.strains[strain].StopGrowth()
            if self.integrator.y[self.numstrains] < self.EmptySubstrateThreshold:
                self.StopGrowth()
            solout(self.integrator.t,self.integrator.y)
        self.AllowGrowth()
        return np.vstack(localtraj)
    
    # growth only needs to final state
    def GrowthSciPyIntegrator(self,initialcells):
        traj = self.Trajectory(initialcells)
        return traj[-1,:self.numstrains]
    

#class GrowthDynamicsTimeIntegrator(GrowthDynamics):
    ## deprecated
    
    #def __init__(self,numstrains = None, **kwargs):
        #if kwargs.get("mixingtime") is None:
            #kwargs["mixingtime"] = 24.
        #super(GrowthDynamicsTimeIntegrator,self).__init__(numstrains = numstrains,**kwargs)

        #self.dyn = TimeIntegrator(dynamics = self.f,initialconditions = np.ones(self.numstrains + 1), params = None)
        #self.dyn.SetEndCondition("maxtime",self.env.mixingtime)
        
        ##self.dyn.SetEndCondition("reachzero",self.numstrains)

        #self.x = np.zeros(self.numstrains)
        #self.otherinitialconditions = np.array([self.env.substrate])

        
    ## simplest dynamics of pure growth and resource use
    ## implemented already in much faster way in parent class
    #def f(self,t,x,params):
        #return np.concatenate([ self.growthrates * x[:self.numstrains],
                                #np.array([ -np.sum(self.growthrates / self.yieldfactors * x[:self.numstrains]) ])
                              #])
    
    ## base growth function to use for time integrator dynamics
    #def Growth(self,initialcells = None):
        #ic = self.checkInitialCells(initialcells)
        #ic = np.concatenate([ic,self.otherinitialconditions])
        #self.dyn.ResetInitialConditions(ic)
        #final = self.dyn.IntegrateToEndConditions()
        #print ic,final
        #return final[:self.numstrains]


    ## should also work for more complicated dynamics implemented in classes inherited from this one
    #def Trajectory(self,timestep,initialconditions):
        ## helper routine
        #def AddTimeToOutputVector(x,t):
            #return np.array([np.concatenate([np.array([t]),x])],dtype=np.float)
        
        ## set initial conditions
        #initialconditions[:self.numstrains] = self.checkInitialCells(initialconditions[:self.numstrains])
        #if len(initialconditions) >= self.numstrains:
            #initialconditions = np.concatenate([initialconditions,self.otherinitialconditions])
        #self.dyn.ResetInitialConditions(initialconditions)
        
        ## generate first entry in output data
        #r = AddTimeToOutputVector(self.dyn.x,0)
        #while not self.dyn.HasEnded():
            #t = self.dyn.IntegrationStep(timestep)
            #r = np.concatenate([r,AddTimeToOutputVector(self.dyn.x,t)],axis=0)
        
        ## return output
        #return r


class GrowthDynamicsPublicGoods(GrowthDynamicsODE):
    def __init__(self,numstrains = None,**kwargs):
        
        #if kwargs.get("mixingtime") is None:
            #kwargs["mixingtime"] = 12.
        
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

        self.otherinitialconditions = np.array([self.env.substrate,0])
        
        self.__onlypositivecoefficients = kwargs.get("onlypositivecoefficients",True)

        
    
    # dynamics for all strains, then substrate, then public good
    def dynamics(self,t,x):
        # public good can influence growth rates and yield
        a = self.GR(x)
        y = self.YF(x)
        return np.concatenate(  [   a*x[:-2],                                           # growth of strains
                                    np.array([  -np.sum(a/y*x[:-2]),                    # decay of nutrients
                                                np.sum(self.__PGProduction*x[:-2])]) ]) # pg
    
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
    
                            
    



class GrowthDynamicsAntibiotics(GrowthDynamicsODE):
    def __init__(self,**kwargs):
        super(GrowthDynamicsAntibiotics,self).__init__(**kwargs)
        
        self.__params = {   'kappa' :         kwargs.get("kappa",2),
                            'gamma' :         kwargs.get("gamma",2),
                            'PGproduction' :  np.array(kwargs.get("PGproduction",np.zeros(self.numstrains)),dtype=np.float64),
                            'PGreductionAB' : kwargs.get("PGreductionAB",1),
                            'PGconc' :        kwargs.get("PGconc",0),   # initial condition PG concentration
                            'ABconc' :        kwargs.get("ABconc",.5)}  # initial concentration antibiotics measured in zMIC

        assert len(self.__params['PGproduction']) == self.numstrains, "PG production not defined correctly"
        assert sum(self.__params['PGproduction']) > 0, "PG is not produced"

        self.otherinitialconditions = np.array([self.env.substrate,self.__params['PGconc'],self.__params['ABconc']])
    
    
    def beta(self,abconc):
        bk = np.power(abconc,self.__params['kappa'])
        return 1 - (1+self.__params['gamma'])*bk/(bk + self.__params['gamma'])
    
    def growthr(self,substrate,abconc):
        if substrate > 0:
            return self.growthrates * self.beta(abconc)
        else:
            return np.zeros(self.numstrains)
    
    def dynamics(self,t,x):
        a = self.growthr(x[-3],x[-1])
        return np.concatenate([ a*x[:-3],                                                      # growth of strains
                                np.array( [ -np.sum(a/self.yieldfactors*x[:-3]),               # decay of nutrients
                                            np.sum(self.__params['PGproduction']*x[:-3]),      # production of public good
                                            -self.__params['PGreductionAB']*x[-1]*x[-2] ])])   # reduction of antibiotics by public good
    
    
    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsAntibiotics,self).ParameterString() +r
        s += "*** antibiotic parameters ***" +r
        s += "  initial concentration " + str(self.__params['ABconc']) +r
        s += "  gamma                 " + str(self.__params['gamma']) +r
        s += "  kappa                 " + str(self.__params['kappa']) +r
        s += "  enzyme production     " + self.arraystring(self.__params['PGproduction']) +r
        s += "  enzyme activity       " + str(self.__params['PGreductionAB']) +r
        s += "  enzyme initial conc.  " + str(self.__params['PGconc']) +r
        return s



class GrowthDynamicsAntibiotics2(GrowthDynamicsODE):
    def __init__(self,**kwargs):
        super(GrowthDynamicsAntibiotics2,self).__init__(**kwargs)
        
        self.__params = {   'kappa' :         kwargs.get("kappa",2),
                            'gamma' :         kwargs.get("gamma",2),
                            'ProductionEfficiency' :  np.array(kwargs.get("AB_Production_Efficiency",np.zeros(self.numstrains)),dtype=np.float64),
                            'ABconc' :        kwargs.get("ABconc",.5)}  # initial concentration antibiotics measured in zMIC

        assert len(self.__params['ProductionEfficiency']) == self.numstrains, "PG production not defined correctly"
        assert sum(self.__params['ProductionEfficiency']) > 0, "PG is not produced"
        
        self.otherinitialconditions = np.array([self.env.substrate,self.__params['ABconc']])
    
    def beta(self,abconc):
        if abconc >= 1e-10:
            bk = np.power(abconc,self.__params['kappa'])
            return (1. - bk)/(1 + bk/self.__params['gamma'])
        else:
            return 1.
    
    def growthr(self,substrate,abconc):
        # new integration scheme needs a more relaxed version of when substrate is empty
        # this variable is set as a (constant) fraction of the average yield
        if substrate > self.EmptySubstrateThreshold:
            return self.growthrates * self.beta(abconc)
        else:
            return np.zeros(self.numstrains)
    
    def dynamics(self,t,x):
        a  = self.growthr(x[-2],x[-1])
        a0 = self.growthr(x[-2],0)
        return np.concatenate([
                                a*x[:self.numstrains],                                                              # growth of strains
                                np.array([
                                    -np.sum(a0/self.yieldfactors*x[:self.numstrains]),                              # decay of nutrients
                                    -np.sum(self.__params['ProductionEfficiency']*x[:self.numstrains]) * x[-1]      # reduction of antibiotics by cells
                                    ])
                                ])
    
    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsAntibiotics2,self).ParameterString() +r
        s += "*** antibiotic parameters ***" +r
        s += "  initial concentration    " + str(self.__params['ABconc']) +r
        s += "  gamma                    " + str(self.__params['gamma']) +r
        s += "  kappa                    " + str(self.__params['kappa']) +r
        s += "  Enzyme Prod & Efficiency " + self.arraystring(self.__params['ProductionEfficiency']) +r
        return s
        

class GrowthDynamicsAntibiotics3(GrowthDynamicsODE):
    def __init__(self,**kwargs):

        super(GrowthDynamicsAntibiotics3,self).__init__(**kwargs)
        self.__params = {   'kappa' :            kwargs.get("kappa",2),
                            'gamma' :            kwargs.get("gamma",2),
                            'BL_Production':     np.array(kwargs.get("BL_Production",np.zeros(self.numstrains)),dtype=np.float64),
                            'BL_Efficiency':     kwargs.get("BL_Efficiency",1e-2),
                            'AB_Conc' :          kwargs.get("AB_Conc",.5),  # initial concentration antibiotics measured in zMIC
                            'AB_Conc_threshold': kwargs.get("AB_Conc_threshold",1e-10),
                            'AB_Diffusivity':    kwargs.get("AB_Diffusivity",1e-3),
                            'BL_Diffusivity':    kwargs.get("BL_Diffusivity",1e-3),
                            'VolumeSeparation':  kwargs.get("VolumeSeparation",1)
                            }

        assert len(self.__params['BL_Production']) == self.numstrains, "PG production not defined correctly"
        assert sum(self.__params['BL_Production']) > 0, "PG is not produced"
        
        self.otherinitialconditions =   np.concatenate([
                                            np.array([self.env.substrate]),      # substrate
                                            np.zeros(self.numstrains),           # internal enzyme concentrations
                                            np.zeros(self.numstrains),           # internal antibiotics concentrations
                                            np.array([0]),                       # external enzyme concentration
                                            np.array([self.__params['AB_Conc']]) # external antibiotics concentration
                                        ])

    def beta(self,abconc):
        if np.any(abconc >= self.__params['AB_Conc_threshold']):
            bk = np.power(abconc,self.__params['kappa'])
            return (1. - bk)/(1 + bk/self.__params['gamma'])
        else:
            return 1.
    
    def growthr(self,substrate,abconc):
        # new integration scheme needs a more relaxed version of when substrate is empty
        # this variable is set as a (constant) fraction of the average yield
        if substrate > self.EmptySubstrateThreshold:
            return self.growthrates * self.beta(abconc)
        else:
            return np.zeros(self.numstrains)
    
    def dynamics(self,t,x):
        a  = self.growthr(x[self.numstrains],x[2*self.numstrains + 1:3*self.numstrains+1])
        a0 = self.growthr(x[self.numstrains],0)
        return np.concatenate([
                                    a*x[:self.numstrains],                                                              # growth of strains
                                    np.array([-np.sum(a0/self.yieldfactors*x[:self.numstrains])]),                      # decay of nutrients
                                    self.__params['BL_Production'] * x[:self.numstrains] - self.__params['BL_Diffusivity'] * (x[self.numstrains + 1:2 * self.numstrains + 1] - x[-2]),
                                    -self.__params['BL_Efficiency'] * x[self.numstrains + 1:2*self.numstrains + 1] * x[2*self.numstrains + 1:3*self.numstrains + 1] - self.__params['BL_Diffusivity'] * (x[2*self.numstrains + 1:3*self.numstrains + 1] - x[-1]),
                                    np.array([self.__params['BL_Diffusivity'] * np.sum(x[:self.numstrains] * (x[self.numstrains + 1:2*self.numstrains + 1] - x[-2]))]),
                                    np.array([self.__params['AB_Diffusivity'] * np.sum(x[:self.numstrains] * (x[2*self.numstrains + 1:3*self.numstrains + 1] - x[-1])) - self.__params['BL_Efficiency'] * x[-1] * x[-2]])
                                ])

    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsAntibiotics3,self).ParameterString() +r
        s += "*** antibiotic parameters ***" +r
        s += "  Antibiotics Initial Conc    " + str(self.__params['AB_Conc']) +r
        s += "  gamma                       " + str(self.__params['gamma']) +r
        s += "  kappa                       " + str(self.__params['kappa']) +r
        s += "  Enzyme Production           " + self.arraystring(self.__params['BL_Production']) +r
        s += "  Enzyme Efficiency           " + str(self.__params['BL_Efficiency']) +r
        s += "  Enzyme Diffusity rate       " + str(self.__params['BL_Diffusivity']) +r
        s += "  Antibiotic Diffusity rate   " + str(self.__params['AB_Diffusivity']) +r
        s += "  Volume/Timescale Separation " + str(self.__params['VolumeSeparation']) +r
        return s


class GrowthDynamicsAntibiotics4(GrowthDynamics):
    def __init__(self,**kwargs):

        super(GrowthDynamicsAntibiotics4,self).__init__(**kwargs)
        self.__params = {   'kappa' :            kwargs.get("kappa",2),
                            'gamma' :            kwargs.get("gamma",2),
                            'BL_Production':     np.array(kwargs.get("BL_Production",np.zeros(self.numstrains)),dtype=np.float64),
                            'BL_Efficiency':     kwargs.get("BL_Efficiency",1e-2),
                            'AB_Conc' :          kwargs.get("AB_Conc",.5),  # initial concentration antibiotics measured in zMIC
                            'AB_Conc_threshold': kwargs.get("AB_Conc_threshold",1e-10),
                            'AB_Diffusivity':    kwargs.get("AB_Diffusivity",1e-3),
                            'BL_Diffusivity':    kwargs.get("BL_Diffusivity",1e-3),
                            'VolumeSeparation':  kwargs.get("VolumeSeparation",1e-10)
                            }
        
        
        # commonly used parameters for the dynamics
        self.rhosigmaE = self.__params['BL_Production'] / self.__params['BL_Diffusivity']
        self.epssigmaB = self.__params['BL_Efficiency'] / self.__params['AB_Diffusivity']
        self.etarho    = self.__params['BL_Production'] * self.__params['VolumeSeparation']
        self.etasigmaB = self.__params['AB_Diffusivity'] * self.__params['VolumeSeparation']
        

        assert len(self.__params['BL_Production']) == self.numstrains, "PG production not defined correctly"
        assert sum(self.__params['BL_Production']) > 0, "PG is not produced"
        
        self.otherinitialconditions =   np.concatenate([
                                            np.array([self.env.substrate]),      # substrate
                                            np.array([0]),                       # external enzyme concentration
                                            np.array([self.__params['AB_Conc']]) # external antibiotics concentration
                                        ])
        
    def beta(self,abconc):
        if np.any(abconc >= self.__params['AB_Conc_threshold']):
            bk = np.power(abconc,self.__params['kappa'])
            return (1. - bk)/(1 + bk/self.__params['gamma'])
        else:
            return 1.
    
    def growthr(self,substrate,abconc):
        # new integration scheme needs a more relaxed version of when substrate is empty
        # this variable is set as a (constant) fraction of the average yield
        if substrate > self.EmptySubstrateThreshold:
            return self.growthrates * self.beta(abconc)
        else:
            return np.zeros(self.numstrains)
    
    def dynamics(self,t,x):
        # adabatic approximation (internal concentrations are 'enslaved' to dynamics of outer concentrations)
        enzyme_internal      = x[self.numstrains + 1] + self.rhosigma # external concentration + how much is hold back from production
        antibiotics_internal = x[self.numstrains + 2] / (1. + self.epssigma * (self.rhosigma + enzyme_internal)) # reduction of external concentration due to internal reduction
        
        a   =  self.growthr(x[self.numstrains],antibiotics_internal)
        a0y = -self.growthr(x[self.numstrains],0)/self.yieldfactors
        
        return np.concatenate([
                    a * x[:self.numstrains],    # cell growth
                    np.array([
                            np.sum(a0y * x[:self.numstrains]), # depletion of nutrients
                            np.sum(self.etarho * x[:self.numstrains]), # production of enzyme
                            -self.__params['BL_Efficiency'] * x[self.numstrains + 1] * x[self.numstrains + 2] + self.etasigmaB * np.dot(x[:self.numstrains],antibiotics_internal - x[self.numstrains + 2]) # reduction of antibiotics
                            ])
                    ])


    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsAntibiotics4,self).ParameterString() +r
        s += "*** antibiotic parameters ***" +r
        s += "  Antibiotics Initial Conc    " + str(self.__params['AB_Conc']) +r
        s += "  gamma                       " + str(self.__params['gamma']) +r
        s += "  kappa                       " + str(self.__params['kappa']) +r
        s += "  Enzyme Production           " + self.arraystring(self.__params['BL_Production']) +r
        s += "  Enzyme Efficiency           " + str(self.__params['BL_Efficiency']) +r
        s += "  Enzyme Diffusity rate       " + str(self.__params['BL_Diffusivity']) +r
        s += "  Antibiotic Diffusity rate   " + str(self.__params['AB_Diffusivity']) +r
        s += "  Volume/Timescale Separation " + str(self.__params['VolumeSeparation']) +r
        return s



class GrowthDynamicsPyoverdin(GrowthDynamicsODE):
    def __init__(self,**kwargs):
        super(GrowthDynamicsPyoverdin,self).__init__(**kwargs)
        
        self.__params = {  'PVDproduction' :  np.array(kwargs.get("PVDproduction",np.zeros(self.numstrains)),dtype=np.float64),
                            'PVDincreaseS'  :  kwargs.get("PVDincreaseS",1),
                            'PVDmaxFactorS' :  kwargs.get("PVDmaxFactorS",1),   # initial condition PG concentration
                            'PVDconc' :        kwargs.get("PVDconc",0)}  # initial concentration antibiotics measured in zMIC

        assert len(self.__params['PVDproduction']) == self.numstrains, "PVD production not defined correctly"
        assert sum(self.__params['PVDproduction']) > 0, "PVD is not produced"
        
        self.otherinitialconditions = np.array([self.env.substrate,self.__params['PVDconc'],self.env.substrate])


    def dynamics(self,t,x):
        p = self.__params['PVDincreaseS'] if x[-1] <= self.env.substrate * self.__params['PVDmaxFactorS'] else 0
        if x[-3] >= 0:
            a = self.growthrates
        else:
            a = np.zeros(self.numstrains)
        return np.concatenate([ a*x[:-3],   np.array([  np.sum(-a*x[:-3]/self.yieldfactors) + p * x[-2],
                                                        np.sum(self.__params['PVDproduction']*x[:-3]),
                                                        p * x[-2]   ])])

    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsPyoverdin,self).ParameterString() +r
        s += "*** pyoverdin parameters ***" +r
        s += "  initial concentration " + str(self.__params['PVDconc']) +r
        s += "  increase rate S       " + str(self.__params['PVDincreaseS']) +r
        s += "  max ratio S           " + str(self.__params['PVDmaxFactorS']) +r
        s += "  pyoverdin production  " + self.arraystring(self.__params['PVDproduction']) +r
        return s



class GrowthDynamicsPyoverdin2(GrowthDynamicsODE):
    def __init__(self,**kwargs):

        super(GrowthDynamicsPyoverdin2,self).__init__(**kwargs)
        
        self.__params = {  'PVDproduction' :  np.array(kwargs.get("PVDproduction",np.zeros(self.numstrains)),dtype=np.float64),
                            'PVDincreaseS'  :  kwargs.get("PVDincreaseS",1),
                            'PVDmaxFactorS' :  kwargs.get("PVDmaxFactorS",1),   # initial condition PG concentration
                            'PVDconc' :        kwargs.get("PVDconc",0)}  # initial concentration antibiotics measured in zMIC

        assert len(self.__params['PVDproduction']) == self.numstrains, "PVD production not defined correctly"
        assert sum(self.__params['PVDproduction']) > 0, "PVD is not produced"
        
        self.otherinitialconditions = np.array([self.env.substrate,self.__params['PVDconc'],self.env.substrate])


    def dynamics(self,t,x):
        p = self.__params['PVDincreaseS'] if x[-1] <= self.env.substrate * self.__params['PVDmaxFactorS'] else 0
        if x[-3] > 0:
            a = self.growthrates
        else:
            a = np.zeros(self.numstrains)
        return np.concatenate([ a*x[:-3],   np.array([  np.sum(-a*x[:-3]/self.yieldfactors) + p * x[-2],
                                                        np.sum(self.__params['PVDproduction']*x[:-3]),
                                                        p * x[-2]   ])])

    
    def ParameterString(self):
        r = '\n'
        s  = super(GrowthDynamicsPyoverdin2,self).ParameterString() +r
        return s




class GrowthDynamicsPyoverdin3(GrowthDynamicsODE):
    def __init__(self,**kwargs):
        super(GrowthDynamicsPyoverdin3,self).__init__(**kwargs)
        
        self.__params = {  'InternalIronYieldCoefficient' : np.array(kwargs.get("PVD_Internal_Yield",np.ones(self.numstrains)),dtype=np.float64),
                            'Production' :                   np.array(kwargs.get("PVD_Production",np.zeros(self.numstrains)),dtype=np.float64),
                            'InitialInternalIron' :          np.array(kwargs.get("PVD_Initial_Internal_Iron",np.zeros(self.numstrains)),dtype=np.float64),
                            'MatchingReceptors' :            np.array(kwargs.get("PVD_Matching_Receptors",np.zeros(self.numstrains)),dtype=np.float64),
                            'BaseIronInflux':                         kwargs.get("PVD_Base_Iron_Influx",1),
                            'Kpvd' :                                  kwargs.get("PVD_Kpvd",1e-30),
                            'TotalIron' :                             kwargs.get("PVD_Total_Iron",1e3),
                            'Efficiency' :                            kwargs.get("PVD_Efficiency",1e-3),
                            'YieldDependence' :                       kwargs.get("PVD_Yield_Dependence","linear")
                        }
        
        assert len(self.__params['Production']) == self.numstrains, "PVD production not defined correctly"
        assert np.sum(self.__params['Production']) > 0, "PVD is not produced"
        
        if self.__params['YieldDependence'] == 'linear':
            self.IronYield = self.IronYieldLinear
        elif self.__params['YieldDependence'] == 'exp':
            self.IronYield = self.IronYieldExp
        else:
            raise NotImplementedError
        self.otherinitialconditions = np.concatenate([self.__params['InitialInternalIron'],np.array([self.env.substrate])])
        
    
    
    def g(self,iron,pvd):
        r = 0
        if pvd > 0:
            a = (iron + pvd - 1)/(2.*pvd)
            b = iron/pvd
            if a*a > b:
                r = a - np.sqrt(a*a - b)
        return r
    
    def IronYieldExp(self,x):
        return self.yieldfactors * (1 - np.exp(-self.__params['InternalIronYieldCoefficient'] * x))
    
    def IronYieldLinear(self,x):
        return self.yieldfactors + self.__params['InternalIronYieldCoefficient'] * x
    
    def dynamics(self,t,x):
        y = self.IronYield( x[self.numstrains:2*self.numstrains] )
        totalPVD = np.sum(self.__params['Production']/self.growthrates * x[:self.numstrains])
        totalPopSize = np.sum(x[:self.numstrains])
        if totalPopSize > 0:
            pvdFe = self.g(self.__params['TotalIron']/self.__params['Kpvd'],totalPVD/self.__params['Kpvd']) * totalPVD / totalPopSize
        else:
            pvdFe = 0.
        if x[-1] > 0:
            a,ay = np.transpose(np.array([[gr,gr/y[i]] if y[i] > 0 else [0.,0.] for i,gr in enumerate(self.growthrates)]))
        else:
            a,ay = np.zeros((2,self.numstrains))
        return np.concatenate([
            a*x[:self.numstrains],
            -a*x[self.numstrains:2*self.numstrains] + self.__params['Efficiency'] * pvdFe + self.__params['BaseIronInflux'],
            np.array([-np.sum(ay * x[:-3])])
            ])
    

    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsPyoverdin3,self).ParameterString() +r
        s += "*** Pyoverdin parameters ***" +r
        s += "  Initial internal iron " + self.arraystring(self.__params['InitialInternalIron']) +r
        s += "  Pyoverdin production  " + self.arraystring(self.__params['Production']) +r
        s += "  Yield effect          " + self.arraystring(self.__params['InternalIronYieldCoefficient']) + r
        s += "  Base Iron influx      " + str(self.__params['BaseIronInflux']) +r
        s += "  Kpvd                  " + str(self.__params['Kpvd']) +r
        s += "  TotalIron             " + str(self.__params['TotalIron']) +r
        s += "  Efficiency            " + str(self.__params['Efficiency']) +r
        return s

        
class GrowthDynamicsPyoverdin4(GrowthDynamicsODE):
    def __init__(self,**kwargs):

        super(GrowthDynamicsPyoverdin4,self).__init__(**kwargs)
        
        self.__params = {  'YieldIncreaseFactor' : kwargs.get("PVD_Yield_Increase_Factor",2),
                            'Production' :          np.array(kwargs.get("PVD_Production",np.zeros(self.numstrains)),dtype=np.float64)
                        }
        
        assert len(self.__params['Production']) == self.numstrains, "PVD production not defined correctly"
        assert np.sum(self.__params['Production']) > 0, "PVD is not produced"
        assert self.__params['YieldIncreaseFactor'] > 0, "Effect on yield not properly defined"

        self.otherinitialconditions = np.array([self.env.substrate])

        
    def dynamics(self,t,x):
        n = np.sum(x[:self.numstrains])
        if n>0:
            y = self.yieldfactors * (self.__params['YieldIncreaseFactor'] - (self.__params['YieldIncreaseFactor'] - 1.)*np.exp(-np.dot(self.__params['Production'],x[:self.numstrains])/n))
        else:
            y = self.yieldfactors
        if x[-1] > 0:
            a = self.growthrates
        else:
            a = np.zeros(self.numstrains)
        return np.concatenate([
                                    a * x[:self.numstrains],
                                    np.array([-np.sum(a * x[:self.numstrains]/y)])
                              ])

                            
    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsPyoverdin4,self).ParameterString() +r
        s += "*** Pyoverdin parameters ***" +r
        s += "  Pyoverdin production  " + self.arraystring(self.__params['Production']) +r
        s += "  Yield effect          " + str(self.__params['YieldIncreaseFactor']) + r
        return s
            
        
    
    
class GrowthDynamicsPyoverdin5(GrowthDynamicsODE):
    def __init__(self,**kwargs):
        super(GrowthDynamicsPyoverdin5,self).__init__(**kwargs)
        
        self.__params = {  'YieldIncreaseFactor' : kwargs.get("PVD_Yield_Increase_Factor",2),
                            'Production' :          np.array(kwargs.get("PVD_Production",np.zeros(self.numstrains)),dtype=np.float64)
                        }
        
        assert len(self.__params['Production']) == self.numstrains, "PVD production not defined correctly"
        assert np.sum(self.__params['Production']) > 0, "PVD is not produced"
        assert self.__params['YieldIncreaseFactor'] > 0, "Effect on yield not properly defined"

        self.otherinitialconditions = np.array([self.env.substrate,0])

        
    def dynamics(self,t,x):
        n = np.sum(x[:self.numstrains])
        if n>0:
            y = self.yieldfactors * (self.__params['YieldIncreaseFactor'] - (self.__params['YieldIncreaseFactor'] - 1.)*np.exp(-x[-1]))
        else:
            y = self.yieldfactors
        if x[-2] > 0:
            a = self.growthrates
        else:
            a = np.zeros(self.numstrains)
        return np.concatenate([
                                    a * x[:self.numstrains],
                                    np.array([-np.sum(a * x[:self.numstrains]/y),
                                              np.dot(self.__params['Production'],x[:self.numstrains]) ])
                              ])

                            
    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsPyoverdin5,self).ParameterString() +r
        s += "*** Pyoverdin parameters ***" +r
        s += "  Pyoverdin production  " + self.arraystring(self.__params['Production']) +r
        s += "  Yield effect          " + str(self.__params['YieldIncreaseFactor']) + r
        return s
            
        
    
    
    
    
class GrowthDynamicsApprox(GrowthDynamics):
    def __init__(self,**kwargs):
        
        super(GrowthDynamicsApprox,self).__init__(**kwargs)
        
        self.__model = kwargs.get('model','GY')
        self.__modelparameters = np.array(kwargs.get('modelparameters',[]),dtype = np.float)
            
        # rewrite parameters for later use
        self.__a  = np.mean(self.growthrates)
        self.__da = (self.growthrates - self.__a)/self.__a
        self.__y  = np.mean(self.yieldfactors)
        self.__dy = (self.yieldfactors - self.__y)/self.__y
        
        self.__sy = self.env.substrate * self.__y
        
    def Growth(self,initialcells = None):
        ic  = self.checkInitialCells(initialcells) # generate list with same dimensions as number of microbial strains
        n = np.sum(ic) * 1.
        if n > 0:
            x = ic/n
            return n * x * np.power(self.__sy/n * self.correction_term(ic[0],ic[1]),1 + self.__da)
        else:
            return np.zeros(2)
    
    def correction_term(self,m1,m2):
        x = 0
        if m1+m2>0:x=float(m1)/(m1+m2)
        r = 1 + self.__dy[0] * (2.*x-1.)
        if self.__model == 'AB':
            if m1 * self.__modelparameters[0] + m2 * self.__modelparameters[1] < 1:
                r = 0
        elif self.__model == 'PVD':
            r *= self.__modelparameters[0]
        return r

    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsApprox,self).ParameterString() +r
        s += "*** Approximation parameters ***" +r
        s += "  Model             " + str(self.__model) + r
        s += "  Model parameters  " + self.arraystring(self.__modelparameters) +r
        return s
        

class GrowthDynamicsResourceExtraction(GrowthDynamicsODE):
    # implement something similar to
    # Elhanati et al, TPB 2011
    # Behar et al, TPB 2014
    # Behar et al, EPL 2015
    # but exclude death rates, as this is replaced by the finite amount of extractable resources and finite time of a cycle
    
    def __init__(self,**kwargs):
        super(GrowthDynamicsResourceExtraction,self).__init__(**kwargs)
        
        self.__params = dict()
        self.__params['ExtractionMaxRate']     = np.array(kwargs.get('ExtractionMaxRate',     np.zeros(self.numstrains)),dtype=np.float)
        self.__params['ExtractionKm']          = np.array(kwargs.get('ExtractionKm',          100. * np.ones(self.numstrains)), dtype=np.float)
        self.__params['InitiallyExtractedRes'] =          kwargs.get('InitiallyExtractedRes', 0)
    
        assert len(self.__params['ExtractionMaxRate']) == self.numstrains, 'Extraction rates not defined for each strain'
        assert len(self.__params['ExtractionKm'])      == self.numstrains, 'ExtractionKm not defined for each strain'
    
        # extractable resources, extracted resources
        self.otherinitialconditions = np.array([self.env.substrate * self.__params['InitiallyExtractedRes'],self.env.substrate * (1-self.__params['InitiallyExtractedRes'])])
        
    def dynamics(self,t,x):
        # growth rates depend linearly on amount of available nutrients
        if x[self.numstrains] > 0:      a    = self.growthrates * x[self.numstrains]
        else:                           a    = np.zeros(self.numstrains)
        
        # extraction dynamics depends on MM kinetics, if extractable resources available
        if x[self.numstrains+1] > 0:    extr = np.sum(x[:self.numstrains]*self.__params['ExtractionMaxRate']/(x[:self.numstrains] + self.__params['ExtractionKm']))
        else:                           extr = 0
        
        return np.concatenate([
            a * x[:self.numstrains],                                     # growth
            np.array([                                      
                extr - np.dot(a/self.yieldfactors, x[:self.numstrains]), # resources are extracted and used for growth
                -np.sum(extr)                        # extractable resources decay
                ])
            ])
    
    def ParameterString(self):
        r  = '\n'
        s  = super(GrowthDynamicsResourceExtraction,self).ParameterString() +r
        s += "*** Resource Extraction parameters ***" +r
        s += "  Extraction rate   " + self.arraystring(self.__params['ExtractionMaxRate']) +r
        s += "  Extraction Km     " + self.arraystring(self.__params['ExtractionKm']) +r
        s += "  Initially Extracted Resources " +  str(self.__params['InitiallyExtractedRes']) +r
        return s
    
    
