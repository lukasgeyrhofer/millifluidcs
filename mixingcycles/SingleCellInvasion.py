#!/usr/bin/env python

import numpy as np
import argparse
import sys

import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-m","--maxM",type=float,default=1e6)
parser = gc.AddGrowthParameters(parser)
args = parser.parse_args()

g = gc.GrowthDynamics(**vars(args))

initialcells1 = np.arange(100)
step = 1
last=100

while last <= args.maxM:
    growth1,growth2 = g.getGrowthMatrix(size = (initialcells1,np.array([1.,0.])))
    for m,n1,n2 in zip(initialcells1,growth1,growth2):
        print "{:6d} {:15.10f} {:15.10f} {:15.10f}".format(m,n1[0],n2[0],n1[1])
    step *= 10
    first = last
    last *= 10
    initialcells1 = np.arange(first,last,step)
    
