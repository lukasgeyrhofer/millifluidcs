#!/usr/bin/env python
import numpy as np
import argparse
import sys,math

import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.addgrowthparameters(parser)
parser.add_argument("-n","--count",type=int,default=1000)
parser.add_argument("-I","--initialconditions",default=[1.,1.])
args = parser.parse_args()


g = gc.StochasticGrowthDynamics(**vars(args))


for i in range(args.count):
    t,n = g.growth(args.initialconditions)
    print("%5d %.6f %5d %5d"%(i,t,n[0],n[1]))
