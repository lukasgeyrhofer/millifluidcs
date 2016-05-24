#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

from growthclasses import growthdynamics
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
parser.add_argument("-k","--droplets",type=int,default=100)
parser.add_argument("-a","--growthrates",nargs=2,default=np.array([1.,1.]))
parser.add_argument("-y","--yieldrates",nargs=2,default=np.array([.90,1.]))
parser.add_argument("-s","--substrate",default=1e5,type=float)
parser.add_argument("-d","--dilution",type=float,default=2e-4)
parser.add_argument("-T","--mixingtime",type=float,default=24.)

parser.add_argument("-N","--initialsize",nargs=2,default=np.array([1e5,1e5]))
parser.add_argument("-m","--mixingcycles",type=int,default=20)
parser.add_argument("-e","--epsilson",type=int,default=1/60.)

parser.add_argument("-p","--production",type=float,default=1e-1)
parser.add_argument("-H","--halflife",type=float,default=10)
parser.add_argument("-A","--epalpha",type=float,default=1)
parser.add_argument("-M","--epmu",type=float,default=1e3)
parser.add_argument("-S","--epsigma",type=float,default=1e2)

args = parser.parse_args()

# generate parameter dictionary
global p
p = {}
p['yield']      = np.array(args.yieldrates,dtype = float)
p['growth']     = np.array(args.growthrates,dtype = float)
p['substrate']  = args.substrate
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilution
p['production'] = args.production
p['decay']      = np.log(2.)/args.halflife

global epp
epp = {}
epp['alpha'] = args.epalpha
epp['mu']    = args.epmu
epp['sigma'] = args.epsigma



g = growthdynamics(growthrates = np.array(args.growthrates,dtype = float), yieldrates = np.array(args.yieldrates,dtype = float), mixingtime = args.mixingtime, dilution = args.dilution, substrate = args.substrate)

fp0 = g.dilution/(1 - g.dilution)*g.substrate*g.yieldrates[0]
maxn = int(10 * fp0)
tmp,growthplain = g.getGrowthMatrix((np.array([0]),np.arange(maxn)))

#print growthplain

y = np.zeros(4)
growthpyoverdine = np.zeros(maxn)
for n in range(maxn):
    y[0] = float(n)
    y[1] = 0.
    y[2] = 0.
    y[3] = g.substrate
    t = 0
    while t <= p['mixingtime']:
        t += args.epsilson
        y = RungeKutta4(dy,y,t,args.epsilson)
    growthpyoverdine[n] = g.dilution * y[0]

for nmean in np.arange(0,maxn,0.2):
    if nmean > 0:
        p = poisson.pmf(np.arange(maxn),nmean)
    else:
        p = np.zeros(maxn)
        p[0] = 1
    print nmean,np.dot(p,growthplain[0]),np.dot(p,growthpyoverdine)





