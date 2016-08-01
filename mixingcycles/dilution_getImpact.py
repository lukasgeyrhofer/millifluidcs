#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.addgrowthparameters(parser)

parser.add_argument("-n","--fixedpoints",nargs="*")
parser.add_argument("-M","--maxN",type=int,default=10)

args = parser.parse_args()


g = gc.growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)
gm1,gm2 = g.getGrowthMatrix(size = args.maxN)

m = np.arange(args.maxN)
fp = np.array(args.fixedpoints,dtype=float)

px,py = gc.PoissonSeedingVectors(m,fp)

#print gm1

for i in m:
    for j in m:
        p  = px[i]*py[j]
        n1 = px[i]*py[j]*gm1[i,j]
        n2 = px[i]*py[j]*gm2[i,j]
        print "{:4d} {:4d} {:.10f} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(i,j,p,n1,n2,fp[0],fp[1])
    print








