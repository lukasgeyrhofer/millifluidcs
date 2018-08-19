#!/usr/bin/env python 

import argparse
import sys,math
import numpy as np

sys.path.append(sys.path[0] + "/..")
import growthclasses as gc


def branches(n,a,y,t,s,d):
    r = d*(s*y+n)
    r[n < s*y/(np.exp(a*t)-1)] = n[n < s*y/(np.exp(a*t)-1)]*np.exp(a*t)*d
    return r


parser = argparse.ArgumentParser()
parser.add_argument("-a","--growthrate",type=float,default=1)
parser.add_argument("-y","--yieldfactor",type=float,default=1)
parser.add_argument("-S","--substrate",type=float,default=1e4)
parser.add_argument("-D","--dilutionrate",type=float,default=2e-4)
parser.add_argument("-t","--TmixMin",type=float,default=5)
parser.add_argument("-T","--TmixMax",type=float,default=20)
parser.add_argument("-d","--dTmix",type=float,default=.1)
parser.add_argument("-m","--maxM",type=int,default=100)
parser.add_argument("-n","--maxN",type=float,default=3)
parser.add_argument("-N","--dN",type=float,default=.1)
args = parser.parse_args()



m = np.arange(args.maxM)
n = np.arange(0,args.maxN+args.dN,args.dN)
tmix = np.arange(args.TmixMin,args.TmixMax+args.dTmix,args.dTmix)

params = {  'growthrates':            np.array([args.growthrate]),
            'yieldfactors':           np.array([args.yieldfactor]),
            'substrateconcentration': args.substrate,
            'dilution':               args.dilutionrate,
            'mixingtime':             args.TmixMin}

g = gc.GrowthDynamics(**params)

p = gc.PoissonSeedingVectors(m,n)

for i in range(len(tmix)):
    g.setMixingTime(tmix[i])
    growth = g.getGrowthVector(args.maxM)
    nn = np.dot(p,growth)
    
    for a,b,c in zip(n,nn,branches(n,args.growthrate,args.yieldfactor,tmix[i],args.substrate,args.dilutionrate)):
        print a,b,c,tmix[i]
    print
    
    




