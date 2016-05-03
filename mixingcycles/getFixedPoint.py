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
        
    mn = np.zeros((size,size,2))
    for i in range(size):
        for j in range(size):
            t = get_time_to_substrate_depletion(i,j)
            mn[i,j] = np.array([i,j]) * np.exp(p['growth']*t) * p['dilution']

    return mn

def prob(m,n):
    if n[0] > 0:
        px = poisson.pmf(m,n[0])
        px[-1] += (1. - np.sum(px))
    else:
        px = np.zeros(len(m))
        px[0] = 1
    if n[1] > 0:
        py = poisson.pmf(m,n[1])
        py[-1] += (1. - np.sum(py))
    else:
        py = np.zeros(len(m))
        py[0] = 1
    return px,py


parser = argparse.ArgumentParser()
parser.add_argument("-m","--maxM",type=int,default=300)
parser.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
parser.add_argument("-Y","--yieldfactor",type=float,nargs="*",default=[1.,2.])
parser.add_argument("-S","--substrateconcentration",type=float,default=1e4)
parser.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
parser.add_argument("-T","--mixingtime",type=float,default=12.)
parser.add_argument("-M","--maxiterations",type=int,default=1000)
parser.add_argument("-A","--alpha",type=float,default=1.)
parser.add_argument("-p","--precision",type=float,default=1e-20)
args = parser.parse_args()



# generate parameter dictionary
p = {}
p['yield']      = np.array(args.yieldfactor)
p['growth']     = np.array(args.growthrates)
p['substrate']  = args.substrateconcentration
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilutionfactor

growthm = getGrowthMatrix(p,size = args.maxM)

m = np.arange(args.maxM)
n = p['dilution']/(1-p['dilution'])*p['substrate']*p['yield']

for i in range(args.maxiterations):
    print "{:4d} {:12.8f} {:12.8f}".format(i,n[0],n[1])
    px,py = prob(m,n)
    # construct iteration function
    fn = np.dot(px,np.dot(py,growthm))
    if np.sum((n-fn)**2) < args.precision:
        break
    n = fn
    

