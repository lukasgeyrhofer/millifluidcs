#!/usr/bin/env python3

import argparse
import numpy as np
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,dilution=False)

parser_dilution = parser.add_argument_group(description = "==== Dilution values ====")
parser_dilution.add_argument("-d","--dilutionmin",type=float,default=1e-6)
parser_dilution.add_argument("-D","--dilutionmax",type=float,default=None)
parser_dilution.add_argument("-K","--dilutionsteps",type=int,default=10)
parser_dilution.add_argument("-L","--dilutionlogscale",default = False, action = "store_true")

parser.add_argument("-m","--maxN",type=int,default=50)
parser.add_argument("-v","--verbose",action="store_true",default=False)
args   = parser.parse_args()

if args.dilutionmax is None:
    dlist = np.array([args.dilutionmin])
else:
    if args.dilutionlogscale:
        dlist = np.power(10,np.linspace(start = np.log10(args.dilutionmin),stop = np.log10(args.dilutionmax), num = args.dilutionsteps))
    else:
        dlist = np.linspace(start = args.dilutionmin,stop = args.dilutionmax,num = args.dilutionsteps)


g                    = gc.GrowthDynamics(**vars(args))
m                    = np.arange(args.maxN)
growth_ICm1          = g.getGrowthMatrix(size=(m,np.array([0,1])))
growth_IC1m          = g.getGrowthMatrix(size=(np.array([0,1]),m))

for dilution in dlist:
    g.env.dilution = dilution

    fp                   = g.getSingleStrainFixedPointsPoissonSeeding(size = args.maxN)
    px                   = gc.PoissonSeedingVectors(m,fp,diff=False)

    Egrowth1_ICm1        = np.dot(growth_ICm1[:,1,0],px[0])
    Egrowth1_IC1m        = np.dot(growth_IC1m[1,:,0],px[0])
    Egrowth2_ICm1        = np.dot(growth_ICm1[:,1,1],px[1])
    Egrowth2_IC1m        = np.dot(growth_IC1m[1,:,1],px[1])

    if args.verbose:
        invasiongrowth1  = g.Growth([fp[0],1])
        invasiongrowth2  = g.Growth([1,fp[1]])

        if fp[0] > 0:                 
            gamma1inv1   = invasiongrowth1[0]/fp[0]
            gamma1inv2   = invasiongrowth2[0]/fp[0]
        else:
            gamma1inv1   = 0
            gamma1inv2   = 0
        if fp[1] > 0:
            gamma2inv1   = invasiongrowth1[1]/fp[1]
            gamma2inv2   = invasiongrowth2[1]/fp[1]
        else:
            gamma2inv1   = 0
            gamma2inv2   = 0

        print("{:e} {:f} {:f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(dilution,g.growthrates[1]/g.growthrates[0],g.yieldfactors[1]/g.yieldfactors[0],fp[0],fp[1],Egrowth1_ICm1,Egrowth2_ICm1,Egrowth1_IC1m,Egrowth2_IC1m,invasiongrowth1[0],invasiongrowth1[1],invasiongrowth2[0],invasiongrowth2[1],gamma1inv1,gamma2inv1,gamma1inv2,gamma2inv2))
        #                                                                                                                               1        2                                 3                                   4     5     6             7             8             9             10                 11                 12                 13                 14         15         16         17
        # invasibility                                                                                                                                                                                                                                                                                        ^strain 1 invaded  ^strain 2 invaded

    else:
        # reduced output. print only fixed points and their relevant stability component
        print("{:e} {:f} {:f} {:.6f} {:6f} {:.6f} {:.6f}".format(dilution,g.growthrates[1]/g.growthrates[0],g.yieldfactors[1]/g.yieldfactors[0],fp[0],fp[1],Egrowth2_ICm1,Egrowth1_IC1m))


