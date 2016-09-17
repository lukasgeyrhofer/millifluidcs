#!/usr/bin/env python3

import argparse
import numpy as np
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,dilution=True)
parser.add_argument("-m","--maxN",type=int,default=200)
args = parser.parse_args()


g = gc.GrowthDynamics(**vars(args))

m = np.arange(args.maxN)

growth1,growth2 = g.getGrowthMatrix(size=(m,np.array([0,1])))

print(g.getSingleStrainFixedPointsPoissonSeeding(size=args.maxN))



for i in range(args.maxN):
    px = gc.PoissonSeedingVectors(m,np.array([i]))[0]
    print("{:5d} {:.10f} {:.10f}".format(i,growth1[i,1],np.dot(px,growth1[:,1])))


