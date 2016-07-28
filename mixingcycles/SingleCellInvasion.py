#!/usr/bin/env python

import numpy as np
import argparse
import sys

from growthclasses import growthdynamics
from growthclasses import addgrowthparamters


parser = argparse.ArgumentParser()
parser.add_argument("-m","--maxM",type=float,default=1e6)
parser = addgrowthparamters(parser)
args = parser.parse_args()

g     = growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)

m = np.arange(1000)
step = 1
last=1000
while last <= args.maxM:
    n0,n1 = g.getGrowthMatrix(size = (m,np.array([1.,0.])))
    for data in zip(m,n0,n1):
        print "{:6d} {:15.10f} {:15.10f} {:15.10f}".format(data[0],data[1][0],data[2][0],data[1][1])
    step *= 10
    first = last
    last *= 10
    m = np.arange(first,last,step)
    
