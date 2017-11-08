#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
sys.path.append(sys.path[0] + '/../')
import growthclasses as gc

from scipy.special import lambertw

def lw(x,onlyReal = True):
    if onlyReal:
        return float(np.real(lambertw(x)))
    else:
        return lambertw(x)


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser.add_argument("-d","--dilutionmin",  type=float, default = 4e-6)
parser.add_argument("-D","--dilutionmax",  type=float, default = 1e-4)
parser.add_argument("-K","--dilutionstep", type=float, default = 2e-6)
parser.add_argument("-L","--dilutionLOG",              default = False, action = "store_true")

parser.add_argument("-m","--maxM",         type=int,   default = 50)
parser.add_argument("-s","--slopeoffset",  type=float, default = .1)
parser.add_argument("-v","--verbose",                  default = False, action = "store_true")

args = parser.parse_args()

g     = gc.GrowthDynamics(**vars(args))
gm    = g.getGrowthMatrix(args.maxM)
m     = np.arange(args.maxM)

if args.dilutionLOG:
    dlist = np.power(10,np.arange(start = args.dilutionmin,stop = args.dilutionmax, step = args.dilutionstep))
else:
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
    nullcline_intersection1 = np.interp(0,(nx_onestep - nx)[::-1],nx[::-1]) # list of values needs to be increasing for numpy to work.
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
    slope1_pullexpectation = (1 - g10[0] + fp[0])/(g11[0] - fp[0])
    slope2_pullexpectation = (g21[1] - fp[1])/(1 - g20[1] + fp[1])
    
    # gamma approximations
    a             = g.growthrates[1]/g.growthrates[0]
    sy            = g.env.substrate * g.yieldfactors
    y             = g.yieldfactors[1]/g.yieldfactors[0]
    fp_appr       = dilution/(1-dilution) * sy
    gamma1_fp_1   = 1 - 1./sy[1] * (np.power(sy[0]/fp_appr[0] +1. - (np.power(sy[0]/fp_appr[0]+1,a)-1)/(fp_appr[0] * y),a)-1.)
    gamma2_1_fp   = fp_appr[1]/sy[1] * (np.power(sy[0]+1.+(1.-np.power(fp_appr[1],1./(1.-a)))/y,a)-1.)
    slope1_approx = 1./((gamma1_fp_1 - 1.)*fp_appr[0])
    slope2_approx = (gamma2_1_fp - 1.)* fp_appr[1]
    
    # new approximations from LW-func
    LWgamma1_fp_1_O1   = 1 - np.power(sy[0]/fp_appr[0] + 1,a)/sy[1]
    LWgamma1_fp_1_O2   = 1 - np.power(sy[0]/fp_appr[0] + 1,a)/sy[1] * (1 + (.5-a)*np.power(sy[0]/fp_appr[0]+1,a)/sy[1])
    tmp_complex        = np.exp(lambertw(-(1-a)*np.power(sy[0]/fp_appr[0],a)/sy[1], k = 0)/(1-a))
    LWgamma1_fp_1_Oe   = float(np.real(tmp_complex))
    LWgamma1_fp_1_Oe_i = float(np.imag(tmp_complex))
    
    slope1_approxLW1 = 1./((LWgamma1_fp_1_O1-1)*fp_appr[0])
    slope1_approxLW2 = 1./((LWgamma1_fp_1_O2-1)*fp_appr[0])
    slope1_approxLWe = 1./((LWgamma1_fp_1_Oe-1)*fp_appr[0])
    
    
    # even newer approximations using LW
    # invasion of strain 2 at SSFP1
    E1   = g.env.substrate * g.yieldfactors[0] / fp_appr[0]
    E2   = g.env.substrate * g.yieldfactors[1] / 1.
    iam1 = 1./(a-1.)
    gamma2_NA = 1. - iam1 * E1/E2 / lw(E2* iam1 * np.exp((1+E2/E1)*iam1))
    # invasion of strain 1 at SSFP2
    E1   = g.env.substrate * g.yieldfactors[0] / 1.
    E2   = g.env.substrate * g.yieldfactors[1] / fp_appr[1]
    iam1 = 1./(a-1.)
    gamma1_NA = iam1 * E1/E2 / lw(E2* iam1 * np.exp((1+E2/E1)*iam1))
    
    
    # exact invasion growth
    invg_fp1 = np.dot(gm[:,1,1],px[0])*dilution
    invg_fp2 = np.dot(gm[1,:,0],px[1])*dilution
    
    

    if args.verbose:
        print '{:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(dilution,-slope1,-slope2,-realslope1,-realslope2,nullcline_intersection1,nullcline_intersection2,fp[0],fp[1],inv12,inv21,excessgrowth1,excessgrowth2,-slope1_pullexpectation,-slope2_pullexpectation,-slope1_approx,-slope2_approx)
    else:
        print '{:14.6e}'.format(dilution),
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(-slope1,-slope2,-slope1_numerics,-slope2_numerics),
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(-slope1_pullexpectation,-slope2_pullexpectation,-slope1_approx,-slope2_approx),
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(-slope1_approxLW1,-slope1_approxLW2,-slope1_approxLWe,LWgamma1_fp_1_Oe_i),
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(gamma1_NA,gamma2_NA,gamma1_NA * fp_appr[0],gamma2_NA * fp_appr[1]),
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(g11[1],g21[0],invg_fp1,invg_fp2)
    
