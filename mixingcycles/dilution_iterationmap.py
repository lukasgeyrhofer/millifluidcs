#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.misc import factorial
from scipy.stats import poisson

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,dilution=True)

parser.add_argument("-m","--maxsize",type=int,default=100)
parser.add_argument("-n","--outputmax",type=float,default=10)
parser.add_argument("-d","--outputdx",type=float,default=.1)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-P","--poissonseeding",default=False,action="store_true")

args = parser.parse_args()

g = gc.GrowthDynamics(**vars(args))
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
            px,py = gc.PoissonSeedingVectors(m,np.array((x,y)))
            gx = np.dot(py,np.dot(px,gm1))
            gy = np.dot(py,np.dot(px,gm2))
        else:
            gx = gm1[x,y]
            gy = gm2[x,y]
        print >> fp,"{} {} {} {}".format(x,y,gx,gy)
    print >> fp
