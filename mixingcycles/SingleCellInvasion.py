#!/usr/bin/env python

import numpy as np
import argparse
import sys

from growthclasses import growthdynamics

parser = argparse.ArgumentParser()
parser.add_argument("-m","--minM",type=int,default=1)
parser.add_argument("-M","--maxM",type=int,default=100)
parser.add_argument("-D","--stepM",type=int,default=1)

parser.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
parser.add_argument("-Y","--yieldrates",type=float,nargs="*",default=[1.,2.])
parser.add_argument("-S","--substrateconcentration",type=float,default=1e4)
parser.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
parser.add_argument("-T","--mixingtime",type=float,default=12.)
args = parser.parse_args()

g     = growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)
m     = np.arange(args.minM,args.maxM,args.stepM)
n0,n1 = g.getGrowthMatrix(size = (m,np.array([1.,0.])))

for data in zip(m,n0,n1):
    print "{:6d} {:14.10f} {:14.10f} {:14.10f}".format(data[0],data[1][0],data[2][0],data[1][1])

