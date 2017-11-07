#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser.add_argument("-n","--initialconditions",nargs="*",type=float)
args = parser.parse_args()


g = gc.GrowthDynamics(**vars(args))

for m1 in args.initialconditions:
    for m2 in args.initialconditions:
        m = np.array([m1,m2])
        print "{:7.2f} {:7.2f} {:7.2f} {:.6e} {:.6e}".format(m1,m2,g.getTimeToDepletion(m),*g.Growth(m))
    print



