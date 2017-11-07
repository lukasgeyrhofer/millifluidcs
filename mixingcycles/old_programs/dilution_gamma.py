#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)
parser.add_argument("-m","--maxN",type=int,default=10)
args = parser.parse_args()


g = gc.GrowthDynamics(**vars(args))

gm1,gm2 = g.getGrowthMatrix(size = args.maxN)

sgm1 = np.repeat(np.transpose(np.array([gm1[:,0]])),args.maxN,axis=1)
sgm2 = np.repeat(np.array([gm2[0,:]]),args.maxN,axis=0)

gamma1 = gm1/sgm1
gamma2 = gm2/sgm2

gamma1[gm1 == 0] = 0
gamma2[gm2 == 0] = 0

for i in range(args.maxN):
    for j in range(args.maxN):
        print("{:4d} {:4d} {:.10f} {:.10f}".format(i,j,gamma1[i,j],gamma2[i,j]))
    print()
