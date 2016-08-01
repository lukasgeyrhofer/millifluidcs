#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.misc import factorial
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

parser.add_argument("-m","--maxsize",type=int,default=300)
parser.add_argument("-n","--outputmax",type=float,default=None)
parser.add_argument("-D","--outputdx",type=float,default=1.)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-P","--poissonseeding",default=False,action="store_true")

args = parser.parse_args()

g = growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)
gm1,gm2 = g.getGrowthMatrix(size = args.maxsize)

m = np.arange(args.maxsize)
if args.poissonseeding:
    outpoints = np.arange(0,args.outputmax,args.outputdx,dtype=float)
else:
    outpoints = np.arange(args.maxsize,dtype=int)
    

writetofile = False
if not args.outfile is None:
    try:
        fp = open(args.outfile,"w")
        writetofile = True
    except:
        fp = sys.stdout
else:
    fp = sys.stdout
    

for x in outpoints:
    for y in outpoints:
        if args.poissonseeding:
            px,py = prob(m,np.array((x,y)))
            gx = np.dot(py,np.dot(px,gm1))
            gy = np.dot(py,np.dot(px,gm2))
        else:
            gx = gm1[x,y]
            gy = gm2[x,y]
        print "{} {} {} {}".format(x,y,gx,gy)
    print
