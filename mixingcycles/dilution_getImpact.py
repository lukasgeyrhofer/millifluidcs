#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

from growthclasses import growthdynamics
from growthclasses import addgrowthparameters

def prob(m,n,cutoff = 1e-100):
    if n[0] > 0:
        px = poisson.pmf(m,n[0])
        px[px<cutoff] = 0.
        px[-1] += (1. - np.sum(px))
    else:
        px = np.zeros(len(m))
        px[0] = 1
    if n[1] > 0:
        py = poisson.pmf(m,n[1])
        py[py<cutoff] = 0.
        py[-1] += (1. - np.sum(py))
    else:
        py = np.zeros(len(m))
        py[0] = 1
    return px,py


parser = argparse.ArgumentParser()
parser = addgrowthparameters(parser)

parser.add_argument("-n","--fixedpoints",nargs="*")
parser.add_argument("-M","--maxN",type=int,default=10)

args = parser.parse_args()


g = growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)

gm1,gm2 = g.getGrowthMatrix(size = args.maxN)

m = np.arange(args.maxN)
fp = np.array(args.fixedpoints,dtype=float)

px,py = prob(m,fp)

#print gm1

for i in m:
    for j in m:
        p  = px[i]*py[j]
        n1 = px[i]*py[j]*gm1[i,j]
        n2 = px[i]*py[j]*gm2[i,j]
        print "{:4d} {:4d} {:.10f} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(i,j,p,n1,n2,fp[0],fp[1])
    print








