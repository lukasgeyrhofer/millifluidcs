#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

from growthclasses import RungeKutta4


def increasegrowth(pyoverdine):
    return 1 + epp['alpha'] / (1 + np.exp(-(pyoverdine - epp['mu'])/ epp['sigma']))


def dy(t,y):
    if y[3] <= 0:
        growth = 0.
        y[3] = 0
    else:
        growth = 1.
    return np.array([
        growth * p['growth'][0] * y[0] * increasegrowth(y[2]),\
        growth * p['growth'][1] * y[1] * increasegrowth(y[2]),\
        growth * p['production'] * y[0] - p['decay'] * y[2],\
        -growth * increasegrowth(y[2]) * (p['growth'][0]/p['yield'][0]*y[0] + p['growth'][1]/p['yield'][1]*y[1]) \
    ])



parser = argparse.ArgumentParser()
parser.add_argument("-a","--growthrates",nargs=2,default=np.array([1.,1.]))
parser.add_argument("-Y","--yieldrates",nargs=2,default=np.array([.90,1.]))
parser.add_argument("-s","--substrateconcentration",default=1e5,type=float)
parser.add_argument("-d","--dilutionrate",type=float,default=2e-4)
parser.add_argument("-T","--mixingtime",type=float,default=24.)
parser.add_argument("-N","--initialsize",nargs=2,default=np.array([1e5,1e5]))
parser.add_argument("-p","--production",type=float,default=1e-1)
parser.add_argument("-H","--halflife",type=float,default=10)
parser.add_argument("-k","--droplets",type=int,default=100)
parser.add_argument("-m","--mixingsteps",type=int,default=20)
parser.add_argument("-e","--epsilson",type=int,default=1/60.)

parser.add_argument("-A","--epalpha",type=float,default=.5)
parser.add_argument("-M","--epmu",type=float,default=1e3)
parser.add_argument("-S","--epsigma",type=float,default=1e3)

args = parser.parse_args()


# generate parameter dictionary
global p
p = {}
p['yield']      = np.array(args.yieldrates)
p['growth']     = np.array(args.growthrates)
p['substrate']  = args.substrateconcentration
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilutionrate
p['production'] = args.production
p['decay']      = np.log(2.)/args.halflife

global epp
epp = {}
epp['alpha'] = args.epalpha
epp['mu']    = args.epmu
epp['sigma'] = args.epsigma


dropp = np.ones(args.droplets)*args.initialsize[0]
dropn = np.ones(args.droplets)*args.initialsize[1]
ep    = np.zeros(args.droplets)

for m in range(args.mixingsteps):
    ep[:] = np.mean(ep)*p['dilution']/args.droplets                  # extracellular product to start into new cycle
    poolp = np.sum(dropp)*p['dilution']/args.droplets                # diluted pool producers (poisson parameter!)
    pooln = np.sum(dropn)*p['dilution']/args.droplets                # diluted pool non-producers

    if poolp > 0:   dropp = 1.*poisson.rvs(poolp,size=args.droplets) # reseed with producer cells
    else:           dropp = np.zeros(args.droplets)
    if pooln > 0:   dropn = 1.*poisson.rvs(pooln,size=args.droplets) # reseed with nonproducer cells
    else:           dropn = np.zeros(args.droplets)
    
    for k in range(args.droplets):
        y = np.array([dropp[k],dropn[k],ep[k],p['substrate']])
        t = 0
        while t <= p['mixingtime']:
            t += args.epsilson
            y = RungeKutta4(dy,y,t,args.epsilson)
            
        print "{:5d} {:5d} {:5.0f} {:5.0f} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e}".format(m,k,dropp[k],dropn[k],ep[k],y[0],y[1],y[2],y[3])
        
        dropp[k] = y[0]
        dropn[k] = y[1]
        ep[k]    = y[2]
    
    print
        
