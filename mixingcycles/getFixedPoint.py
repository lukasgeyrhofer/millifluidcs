#!/usr/bin/env python


import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

def getGrowthMatrix(parameters,size = 300):
        def get_time_to_substrate_depletion(m1,m2,alpha=1.,maxsteps = 10000, precision = 1e-10):
        initialcells = np.array([m1,m2])
        # internal function to have all parameters.
        t0 = 0
        if np.sum(initialcells) > 0:
            # initial condition for iteration is assuming only strain with highest expected yield is present
            p = (1.*initialcells/parameters['yield']).argmax() # get ID of strain
            t1 = np.log(parameters['substrate']*parameters['yield'][p]/initialcells[p]+1.)/parameters['growth'][p] # set first estimate of time for single strain
            i = 0
            while ((t1-t0)/t1)**2 > precision**2:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += alpha*(parameters['substrate']-np.sum(initialcells/parameters['yield']*(np.exp(parameters['growth']*t1)-1.)))/(np.sum(initialcells/parameters['yield']*parameters['growth']*np.exp(parameters['growth']*t1)))
                i+=1
                # should not iterate infinitely
                if i > maxsteps:
                    raise ValueError
            return min(t1,parameters['mixingtime'])
        else:
            return 0
        
        m = np.zeros((size,size,2))
        for i in range(size):
            for j in range(size):
                t = get_time_to_substrate_depletion(i,j)
                m[i,j] = np.array([i,j]) * np.exp(p['growthrates']*t) * p['dilution']

        return m

def fn(n):
    px = poisson.pmf(m,n[0])
    py = poisson.pmf(m,n[1])
    ppxy = np.outer(px,py)
    return n - np.sum(ppxy,mnext,axis=(0,1))
    
def Jac(n):
    px = poisson.pmf(m,n[0])
    py = poisson.pmf(m,n[1])
    
    if n[0] > 0:
        pxn = px*(m/n[0]-1)
    else:
        pxn = -px
    if n[1] > 0:
        pyn = py*(m/n[1]-1)
    else:
        pyn = -py

    return np.array([ [1. - np.sum(np.outer(pxn,py)*mnext[:,:,0]),n[0] - np.sum(np.outer(px,pyn)*mnext[:,:,0])],[n[1] - np.sum(np.outer(pxn,py)*mnext[:,:,1]),1. - np.sum(np.outer(px,pyn)*mnext[:,:,1])]])
    

parser = argparse.ArgumentParser()
parser.add_argument("-m","--maxM",type=int,default=300)
parser.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
parser.add_argument("-Y","--yieldfactor",type=float,nargs="*",default=[1.,2.])
parser.add_argument("-S","--substrateconcentration",type=float,default=1e4)
parser.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
parser.add_argument("-T","--mixingtime",type=float,default=12.)
args = parser.parse_args()



# generate parameter dictionary
p = {}
p['yield']      = np.array(args.yieldfactor)
p['growth']     = np.array(args.growthrates)
p['substrate']  = args.substrateconcentration
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilutionfactor

global mnext
mnext = getGrowthMatrix(p,size = args.maxM)
global m
m = np.arange(args.maxM)


n = p['dilution']/(1-p['dilution'])*p['substrate']*p['yield']

for i in range(args.maxiterations):
    print n
    n = n - args.alpha*np.dot(np.inv(Jac(n)),fn(n))

