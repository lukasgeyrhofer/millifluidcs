#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
import growthclasses as gc
from scipy.special import lambertw


def lw(x):
    return np.array(np.real(lambertw(x)),dtype=float)


def gamma1(m1,m2,params):
    if m1 == 0 and m2 == 0:
        return 0
    elif m1 == 0 and m2 > 0:
        return 0
    elif m2 == 0 and m2 > 0:
        return 1
    else:
        E1 = params['SY1']/m1
        E2 = params['SY2']/m2
        iam1 = 1./(params['a'] - 1)
        return iam1 * E1/E2 / lw(iam1 * E2 * np.exp( iam1 * (1 + E2/E1)))


def gamma2(m1,m2,params):
    if m1 == 0 and m2 == 0:
        return 0
    elif m1 == 0 and m2 > 0:
        return 1
    elif m2 == 0 and m2 > 0:
        return 0
    else:
        E2 = params['SY1']/m1
        E1 = params['SY2']/m2
        iam1 = 1./(1./params['a'] - 1)
        return iam1 * E1/E2 / lw(iam1 * E2 * np.exp( iam1 * (1 + E2/E1)))
    


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser.add_argument("-d","--dilutionmin",  type=float, default = 4e-6)
parser.add_argument("-D","--dilutionmax",  type=float, default = 1e-4)
parser.add_argument("-K","--dilutionstep", type=float, default = 2e-6)
parser.add_argument("-L","--dilutionLOG",              default = False, action = "store_true")

parser.add_argument("-m","--maxM",         type=int,   default = 50)
parser.add_argument("-s","--slopeoffset",  type=float, default = .1)
parser.add_argument("-v","--verbose",                  default = False, action = "store_true")

args  = parser.parse_args()

g     = gc.GrowthDynamics(**vars(args))
gm    = g.getGrowthMatrix(args.maxM)
m     = np.arange(args.maxM)

if args.dilutionLOG:
    dlist = np.power(10,np.arange(start = args.dilutionmin,stop = args.dilutionmax, step = args.dilutionstep))
else:
    dlist = np.arange(start = args.dilutionmin,stop = args.dilutionmax, step = args.dilutionstep)

params = {'SY1' : g.env.substrate * g.yieldfactors[0],
          'SY2' : g.env.substrate * g.yieldfactors[1],
          'a'   : g.growthrates[1] / g.growthrates[0]}


for dilution in dlist:
    g.setDilution(dilution)
    fp_exact   = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxM)
    fp_approx1 = g.env.dilution * g.env.substrate * g.yieldfactors / (1 - g.env.dilution)
    fp_approx2 = fp_approx1 - np.exp(-fp_approx1+1) # account for empty droplets
    fp_approx2[fp_approx2 < 0] = 0
    
    px_e  = gc.PoissonSeedingVectors(m,fp_exact)
    px_a1 = gc.PoissonSeedingVectors(m,fp_approx1)
    px_a2 = gc.PoissonSeedingVectors(m,fp_approx2)
    
    G2_ssfp1 = np.dot(gm[:,1,1],px_e[0]) * dilution
    G1_ssfp2 = np.dot(gm[1,:,0],px_e[1]) * dilution
    
    gamma2_ssfp1 = np.dot([1 - gamma1(mm,1,params) for mm in m],px_e[0]) * dilution
    gamma1_ssfp2 = np.dot([gamma1(1,mm,params) for mm in m],px_e[1]) * dilution
    
    G2_ssfp1_approxGamma = gamma2_ssfp1 * fp_exact[1]
    G1_ssfp2_approxGamma = gamma1_ssfp2 * fp_exact[0]
    
    G2_ssfp1_approxInoc  = (1 - gamma1(fp_exact[0],1,params)) * fp_exact[1]
    G1_ssfp2_approxInoc  = (1 - gamma2(1,fp_exact[1],params)) * fp_exact[0]
    
    print '{:14.6e} {:10.6f}'.format(dilution,dilution * params['SY1']),
    print '{:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(fp_exact[0],fp_exact[1],fp_approx1[0],fp_approx1[1],fp_approx2[0],fp_approx2[1]),
    print '{:14.6e} {:14.6e} {:14.6e}'.format(G2_ssfp1,G2_ssfp1_approxGamma,G2_ssfp1_approxInoc),
    print '{:14.6e} {:14.6e} {:14.6e}'.format(G1_ssfp2,G1_ssfp2_approxGamma,G1_ssfp2_approxInoc)
    
    
