#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser.add_argument("-m","--maxM",type=int,default=50)
parser.add_argument("-d","--dilutionmin",type=float,default=2e-6)
parser.add_argument("-D","--dilutionmax",type=float,default=5e-5)
parser.add_argument("-K","--dilutionstep",type=float,default=2e-6)

parser.add_argument("-s","--slopeoffset",type=float,default=.1)

args = parser.parse_args()

g  = gc.GrowthDynamics(**vars(args))
gm = g.getGrowthMatrix(args.maxM)
m  = np.arange(args.maxM)

dlist = np.arange(start = args.dilutionmin,stop = args.dilutionmax, step = args.dilutionstep)

for dilution in dlist:
    g.setDilution(dilution)
    fp = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxM)
    px = gc.PoissonSeedingVectors(np.arange(args.maxM),fp)

    inv21         = np.dot(gm[:,1,0],px[0,:]) * dilution
    inv12         = np.dot(gm[1,:,1],px[1,:]) * dilution
    excessgrowth1 = np.dot(gm[1:,0,0],px[0,:-1]) * dilution
    excessgrowth2 = np.dot(gm[0,1:,1],px[1,:-1]) * dilution

    slope1 = (1.-excessgrowth1 + fp[0])/(inv21 - fp[0])
    slope2 = (inv12-fp[1])/(1.-excessgrowth2 + fp[1])
    
    nx = np.arange(start = 0,stop = fp[0] + 1.,step = args.slopeoffset)
    nx_onestep = np.zeros(len(nx))
    for i in range(len(nx)):
        px = gc.PoissonSeedingVectors(m,[nx[i],args.slopeoffset])
        nx_onestep[i] = np.dot(px[1],np.dot(px[0],gm[:,:,0]))
    
    nullcline1 = np.interp(0,nx,nx_onestep - nx)
    realslope1 = args.slopeoffset/(nullcline1 - fp[0])
    
    ny = np.arange(start = 0,stop = fp[1] + 1.,step = args.slopeoffset)
    ny_onestep = np.zeros(len(ny))
    for i in range(len(ny)):
        py = gc.PoissonSeedingVectors(m,[args.slopeoffset,ny[i]])
        ny_onestep[i] = np.dot(py[1],np.dot(py[0],gm[:,:,1]))

    nullcline2 = np.interp(0,ny,ny_onestep - ny)
    realslope2 = args.slopeoffset/(nullcline2 - fp[1])

    print '{:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(dilution,-slope1,-slope2,fp[0],fp[1],inv12,inv21,excessgrowth1,excessgrowth2,realslope1,realslope2)
    