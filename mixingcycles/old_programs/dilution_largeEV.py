#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser.add_argument("-a","--growthrateratio",default=1,type=float)
parser.add_argument("-D","--dilution",default=2e-4,type=float)
parser.add_argument("-m","--maxM",default=100,type=int)
parser.add_argument("-L","--logSteps",action="store_true",default=False)
parser.add_argument("-S","--substrateconcentration",default=1e4,type=float)
parser.add_argument("-T","--mixingtime",default=1000,type=float)

parser.add_argument("-R","--deltaparameters",default=False,action="store_true")
parser.add_argument("-y","--minyieldratio",default=.4,type=float)
parser.add_argument("-Y","--maxyieldratio",default=4,type=float)
parser.add_argument("-d","--stepyield",default=.01,type=float)


args = parser.parse_args()

m      = np.arange(args.maxM)
params = { "growthrates"            : np.array([1.,args.growthrateratio]),
           "yieldfactors"           : np.ones(2,dtype=np.float),
           "dilution"               : args.dilution,
           "substrateconcentration" : args.substrateconcentration,
           "mixingtime"             : args.mixingtime}


if args.logSteps:
    ylistexp  = np.arange(args.minyieldratio,args.maxyieldratio+args.stepyield,args.stepyield)
    ylist     = np.power(10,ylistexp)
else:
    ylist     = np.arange(args.minyieldratio,args.maxyieldratio+args.stepyield,args.stepyield)

for y in ylist:
    if args.deltaparameters:
        params["yieldfactors"] = np.array([1-y,1+y])
    else:
        params["yieldfactors"][1] = y
    g = gc.GrowthDynamics(**params)
    n = g.getSingleStrainFixedPointsPoissonSeeding(size=args.maxM)

    growth2_atFP1 = g.getGrowthMatrix(size = np.array([m,np.ones(1)]))[:,0,1]
    growth1_atFP2 = g.getGrowthMatrix(size = np.array([np.ones(1),m]))[0,:,0]
    
    p             = gc.PoissonSeedingVectors(m,n)
    stabFP1       = np.dot(p[0],growth2_atFP1)
    stabFP2       = np.dot(p[1],growth1_atFP2)
    
    print("{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(args.growthrateratio,y,stabFP1,stabFP2,n[0],n[1]))

