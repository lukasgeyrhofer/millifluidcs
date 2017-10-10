#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)
parser.add_argument("-m","--maxM",type=int,default=30)

parser.add_argument("-d","--dilutionmin",type=float,default=2e-6)
parser.add_argument("-D","--dilutionmax",type=float,default=5e-5)
parser.add_argument("-K","--dilutionstep",type=float,default=2e-6)

args = parser.parse_args()


g  = gc.GrowthDynamics(**vars(args))
gm = g.getGrowthMatrix(args.maxM)

dlist = np.arange(start = args.dilutionmin,stop = args.dilutionmax + args.dilutionstep, step = args.dilutionstep)

for dilution in dlist:
    g.setDilution(dilution)
    fp = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxM)
    px = gc.PoissonSeedingVectors(np.arange(args.maxM),fp)
    slope1 = (1./dilution - np.dot(gm[1:,0,0] - gm[:-1,0,0],px[0,:-1]))/np.dot(gm[:,1,0] - gm[:,0,0],px[0,:])
    slope2 = np.dot(gm[1,:,1] - gm[0,:,1],px[1,:])/(1./dilution-np.dot(gm[0,1:,1] - gm[0,:-1,1],px[1,:-1]))

    print '{:.6e} {:14.6e} {:14.6e}'.format(dilution,slope1,slope2)
