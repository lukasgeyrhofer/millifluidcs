#!/usr/bin/env python3

import argparse
import sys,math
import numpy as np
import growthclasses as gc

parser = argparse.ArgumentParser()
parser.add_argument("-D","--delta",type=float,default=.2)
parser.add_argument("-M","--maxexp",type=float,default=1)
args = parser.parse_args()


g = gc.GrowthDynamics(dilution = 1.)

for yexp in np.arange(-args.maxexp,args.maxexp+args.delta,args.delta):
    y = np.array([10**yexp,1.])
    g.yieldfactors = y
    for aexp in np.arange(-args.maxexp,args.maxexp+args.delta,args.delta):
        a = np.array([10**aexp,1.])
        g.growthrates = a
        
        n11 = g.getGrowth(np.array([2,0]))
        n22 = g.getGrowth(np.array([0,2]))
        n12 = g.getGrowth(np.array([1,1]))
        
        print("{:10.6f} {:10.6f} {:.12f} {:.12f}".format(y[0],a[0],n12[0]/n11[0],n12[1]/n22[1]))
        
    print()



