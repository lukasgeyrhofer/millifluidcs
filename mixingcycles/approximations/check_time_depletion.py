#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
sys.path.append(sys.path[0] + '/../')
import growthclasses as gc
from scipy.special import lambertw

def lw(x,complex = False):
    r = lambertw(x)
    if complex:
        return r
    else:
        return np.array(r.real,dtype=float)

def ten_to(x):
    return np.power(10.,x)


parser = argparse.ArgumentParser()
parser_growth = parser.add_argument_group(description = "=== Growth parameters ===")
parser_growth.add_argument("-a","--log10minGrowthRate",default=-1,type=float)
parser_growth.add_argument("-A","--log10maxGrowthRate",default=1,type=float)
parser_growth.add_argument("-y","--log10minYield",default=-1,type=float)
parser_growth.add_argument("-Y","--log10maxYield",default=1,type=float)
parser_growth.add_argument("-S","--substrate",default=5e4,type=float)
parser_growth.add_argument("-D","--dilution",default=4e-5,type=float)
parser_growth.add_argument("-k","--log10stepGrowthRate",default=.1,type=float)
parser_growth.add_argument("-K","--log10stepYield",default=.1,type=float)
parser_growth.add_argument("-T","--mixingtime",default=24,type=float)

parser_alg = parser.add_argument_group(description = "=== Algorithm parameters ===")
parser_alg.add_argument("-m","--maxM",default=50,type=int)

args = parser.parse_args()

params = {  'growthrates':  np.array([1,ten_to(args.log10minGrowthRate)]),
            'yieldfactors': np.array([1,ten_to(args.log10minYield)]),
            'substrateconcentration': args.substrate,
            'dilution': args.dilution,
            'mixingtime': args.mixingtime }
            

g = gc.GrowthDynamics(**params)
dSY1 = args.substrate * args.dilution

alist = ten_to(np.arange(start = args.log10minGrowthRate, stop = args.log10maxGrowthRate, step = args.log10stepGrowthRate))
ylist = ten_to(np.arange(start = args.log10minYield,      stop = args.log10maxYield,      step = args.log10stepYield))

for a in alist:
    for y in ylist:
        g.growthrates = np.array([1,a])
        g.yieldfactors = np.array([1,y])
        
        fp = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxM)

        # calculate time to depletion
        tdepl_near1_exact    = g.getTimeToDepletion(initialcells = np.array([fp[0],1]))
        tdepl_near2_exact    = g.getTimeToDepletion(initialcells = np.array([1,fp[1]]))
        
        tdepl_near1_approxFP = g.getTimeToDepletion(initialcells = np.array([dSY1,1]))
        tdepl_near2_approxFP = g.getTimeToDepletion(initialcells = np.array([1,dSY1*y]))
        
        E1_near1             = args.substrate/fp[0]
        E2_near1             = y * args.substrate
        E1_near2             = args.substrate
        E2_near2             = y * args.substrate/fp[1]
        tdepl_near1_approx   = (1+E2_near1/E1_near1)/(1.-a) + lw(E2_near1/(a-1)*np.exp((1+E2_near1/E1_near1)/(a-1)))
        tdepl_near2_approx   = (1+E2_near2/E1_near2)/(1.-a) + lw(E2_near2/(a-1)*np.exp((1+E2_near2/E1_near2)/(a-1)))
        
        # use time to depletion to estimate growth
        px = gc.PoissonSeedingVectors(np.arange(args.maxM),fp)
        GM2_near1 = np.zeros(args.maxM)
        GM1_near2 = np.zeros(args.maxM)
        for i in range(args.maxM):
            GM2_near1[i] = g.Growth(initialcells = np.array([i,1]))[1]
            GM1_near2[i] = g.Growth(initialcells = np.array([1,i]))[0]
        G2_near1 = np.dot(GM2_near1,px[0])
        G1_near2 = np.dot(GM1_near2,px[1])
        
        Gn2_near1             = g.Growth(initialcells = np.array([fp[0],1]))[1]
        Gn1_near2             = g.Growth(initialcells = np.array([1,fp[1]]))[0]
        Gn2_near1_othermethod = np.exp(a * tdepl_near1_exact + np.log(args.dilution))
        Gn1_near2_othermethod = np.exp(    tdepl_near2_exact + np.log(args.dilution))
        
        G2_near1_approx       = np.exp(a * tdepl_near1_approx + np.log(args.dilution))
        G1_near2_approx       = np.exp(    tdepl_near2_approx + np.log(args.dilution))
        
        # output
        print '{:7.4f} {:7.4f}'.format(a,y),
        #                              1 2
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(tdepl_near1_exact, tdepl_near2_exact, tdepl_near1_approxFP, tdepl_near2_approxFP, tdepl_near1_approx, tdepl_near2_approx),
        #                                                                    3                  4                  5                     6                     6                   8
        print '{:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(G2_near1, G1_near2, Gn2_near1, Gn1_near2),
        #                                                  9         10        11         12
        print '{:14.6e} {:14.6e}'.format(G2_near1_approx, G1_near2_approx),
        #                                13               14
        print '{:14.6e} {:14.6e}'.format(Gn2_near1_othermethod, Gn1_near2_othermethod)
        #                                15                     16
    print
    
