#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc

def gamma1(mm1,mm2,params):
    inv = False
    if params.growthrates[1] < params.growthrates[0]: # strain 1 is the fast growing
        m1  = mm1
        m2  = mm2
        a   = params.growthrates[1]/params.growthrates[0]
        y   = params.yieldfactors[1]/params.yieldfactors[0]
        sy1 = params.yieldfactors[0] * params.substrateconcentration
    else:
        inv = True
        m1  = mm2
        m2  = mm1
        a   = params.growthrates[0]/params.growthrates[1]
        y   = params.yieldfactors[0]/params.yieldfactors[1]
        sy1 = params.yieldfactors[1] * params.substrateconcentration
    
    if m1 == 0 or m2 == 0:
        return 1
    else:
        m   = m1/m2
        if m > 1: #m1 >= m2:
            gamma = 1. - m2/(sy1+m1) * (1./y) * (np.power(sy1/m1+1.,a)-1)
        else:
            gamma = 1. - m2/(sy1+m1) * (1./y) * np.power(m,1+1./a)*np.power(sy1/m1 + np.power(m,1./(1-a)),1./a)
        
        if inv:
            return 1. - gamma
        else:
            return gamma

def gamma2(mm1,mm2,params):
    return 1. - gamma1(mm1,mm2,params)


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,dilution=True)
parser.add_argument("-n","--maxGMsize",default=30,type=int)
args   = parser.parse_args()

g   = gc.GrowthDynamics(**vars(args))
gm  = g.getGrowthMatrix(args.maxGMsize)


for x in range(args.maxGMsize):
    for y in range(args.maxGMsize):
        print '{:4d} {:4d} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(x,y,gamma1(x,y,args),gm[x,y,0]/gm[x,0,0],gamma2(x,y,args),gm[x,y,1]/gm[0,y,1])
    print










