#!/usr/bin/env python


import argparse
import numpy as np
import sys,math
import growthclasses as gc


def substr(n,initialcells,gg):
    if n==0:
        return 0
    else:
        return 1/gg.yieldfactors[1]*(np.exp(gg.growthrates[1]*tsat(n,initialcells,gg))-1)


def tsat(n,initialcells,gg):
    if n<1:n=1
    return 1/gg.growthrates[0] * np.log((gg.env.substrate - substr(n-1,initialcells,gg))* gg.yieldfactors[0]/initialcells + 1)


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)
parser.add_argument("-n","--iterationdepth",default=2,type=int)
parser.add_argument("-m","--initialcells",default=10,type=float)
args = parser.parse_args()


g = gc.GrowthDynamics(**vars(args))

g1 = g.Growth([args.initialcells,0])
g2 = g.Growth([args.initialcells,1])

gamma = g2[0]/g1[0]

for n in range(args.iterationdepth):
    itergamma = 1-substr(n,args.initialcells,g)/g.env.substrate
    print("{:4d} {:.10f} {:.10f} {:e}".format(n,itergamma,gamma,1-itergamma/gamma))




