#!/usr/bin/env python
import numpy as np
import argparse
import sys,math

import growthclasses as gc

from scipy.stats import beta,binom

parser = argparse.ArgumentParser()
parser = gc.addgrowthparameters(parser)
parser.add_argument("-n","--count",type=int,default=1000)
parser.add_argument("-I","--initialconditions",default=[1.,1.])
args = parser.parse_args()


g = gc.StochasticGrowthDynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)


for i in range(args.count):
    t,n = g.growth(args.initialconditions)
    print("%5d %.6f %5d %5d"%(i,t,n[0],n[1]))
