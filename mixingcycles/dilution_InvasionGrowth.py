#!/usr/bin/env python3

import argparse
import numpy as np
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,dilution=True)
parser.add_argument("-m","--maxN",type=int,default=10)
args = parser.parse_args()

g         = gc.GrowthDynamics(**vars(args))
m         = np.arange(args.maxN)
fp        = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxN)
g1m1,g1m2 = g.getGrowthMatrix(size=(m,np.array([0,1])))
g2m1,g2m2 = g.getGrowthMatrix(size=(np.array([0,1]),m))

px,dpx = gc.PoissonSeedingVectors(m,fp,diff=True)

Egm1 = np.dot(g1m1[:,1],px[0])
Eg1m = np.dot(g2m2[1,:],px[1])

invasion1 = g.Growth([fp[0],1])
invasion2 = g.Growth([1,fp[1]])

gamma1inv1 = invasion1[0]/fp[0]
gamma1inv2 = invasion2[0]/fp[0]

gamma2inv1 = invasion1[1]/fp[1]
gamma2inv2 = invasion2[1]/fp[1]


print("{:f} {:f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(g.growthrates[1]/g.growthrates[0],g.yieldfactors[1]/g.yieldfactors[0],Egm1,Eg1m,invasion1[0],invasion1[1],invasion2[0],invasion2[1],gamma1inv1,gamma1inv2,gamma2inv1,gamma2inv2))

