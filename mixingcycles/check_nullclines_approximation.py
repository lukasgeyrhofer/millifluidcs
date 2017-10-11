#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser.add_argument("-d","--dilutionmin",  type=float, default = 4e-6)
parser.add_argument("-D","--dilutionmax",  type=float, default = 1e-4)
parser.add_argument("-K","--dilutionstep", type=float, default = 2e-6)

parser.add_argument("-m","--maxM",         type=int,   default = 50)
parser.add_argument("-s","--slopeoffset",  type=float, default = .1)
parser.add_argument("-v","--verbose",                  default = False, action = "store_true")

args = parser.parse_args()

g     = gc.GrowthDynamics(**vars(args))
gm    = g.getGrowthMatrix(args.maxM)

m     = np.arange(args.maxM)
dlist = np.arange(start = args.dilutionmin,stop = args.dilutionmax, step = args.dilutionstep)

for dilution in dlist:
    g.setDilution(dilution)
    fp = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxM)
    px = gc.PoissonSeedingVectors(np.arange(args.maxM),fp)

    # approximation
    inv21         = np.dot(gm[:,1,0],px[0,:]) * dilution
    inv12         = np.dot(gm[1,:,1],px[1,:]) * dilution
    excessgrowth1 = np.dot(gm[1:,0,0],px[0,:-1]) * dilution
    excessgrowth2 = np.dot(gm[0,1:,1],px[1,:-1]) * dilution
    slope1        = (1.-excessgrowth1 + fp[0])/(inv21 - fp[0])
    slope2        = (inv12-fp[1])/(1.-excessgrowth2 + fp[1])

    
    # compute real slopes from iteration map
    # slope of nullcline 1 
    nx         = np.arange(start = 0,stop = fp[0] + 1.,step = args.slopeoffset)
    nx_onestep = np.zeros(len(nx))
    for i in range(len(nx)):
        px            = gc.PoissonSeedingVectors(m,[nx[i],args.slopeoffset])
        nx_onestep[i] = np.dot(px[1],np.dot(px[0],gm[:,:,0]))*dilution
    nullcline_intersection1 = np.interp(0,(nx_onestep - nx)[::-1],nx[::-1])
    slope1_numerics         = args.slopeoffset/(nullcline_intersection1 - fp[0])
    
    # slope of nullcline 2
    ny         = np.arange(start = 0,stop = fp[1] + 1.,step = args.slopeoffset)
    ny_onestep = np.zeros(len(ny))
    for i in range(len(ny)):
        py            = gc.PoissonSeedingVectors(m,[args.slopeoffset,ny[i]])
        ny_onestep[i] = np.dot(py[1],np.dot(py[0],gm[:,:,1]))*dilution
    nullcline_intersection2 = np.interp(0,(ny_onestep - ny)[::-1],ny[::-1])
    slope2_numerics         = (nullcline_intersection2 - fp[1])/args.slopeoffset
    

    # slopes with <G(m)> = G(<m>)
    g11                    = g.Growth(initialcells = np.array([fp[0],1]))
    g10                    = g.Growth(initialcells = np.array([fp[0] + 1,0]))
    g21                    = g.Growth(initialcells = np.array([1,fp[1]]))
    g20                    = g.Growth(initialcells = np.array([0,fp[1]+1]))
    slope1_pullexpectation = (1-g11[0]-fp[0])/(g10[0] - fp[0])
    slope2_pullexpectation = (g20[1] - fp[1])/(1-g21[1]-fp[1])
    
    # gamma approximations
    a             = g.growthrates[1]/g.growthrates[0]
    sy            = g.env.substrate * g.yieldfactors
    y             = g.yieldfactors[1]/g.yieldfactors[0]
    fp_appr       = dilution/(1-dilution) * sy
    slope1_approx = - sy[1] * (1-dilution) / ((np.power(sy[0]/fp_appr[0] + 1.,a)-1)*fp_appr[0])
    slope2_approx = - (1-fp_appr[1]/sy[1] * (np.power(sy[0]+1+(1-np.power(fp_appr[1],1/(1-a)))/y,a)-1))*fp_appr[1]
    

    if args.verbose:
        print '{:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(dilution,-slope1,-slope2,-realslope1,-realslope2,nullcline_intersection1,nullcline_intersection2,fp[0],fp[1],inv12,inv21,excessgrowth1,excessgrowth2,-slope1_pullexpectation,-slope2_pullexpectation,-slope1_approx,-slope2_approx)
    else:
        print '{:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} '.format(dilution,-slope1,-slope2,-slope1_numerics,-slope2_numerics),
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(-slope1_pullexpectation,-slope2_pullexpectation,-slope1_approx,-slope2_approx)
    
