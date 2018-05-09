#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle
import itertools

import growthclasses as gc

def cont(dx,maxdx,step,maxstep):
    ret = True
    if not maxstep is None:
        if step >= maxstep:
            ret = False
    if maxdx*maxdx < dx*dx:
        ret = False
    return ret


parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-i","--infile",required=True)
parser_io.add_argument("-v","--verbose",action="store_true",default=False)

parser_m = parser.add_argument_group("==== Minimizer Parameters ====")
parser_m.add_argument("-m","--maxsteps",default=None,type=int)
parser_m.add_argument("-p","--precision",default=1e-10,type=float)
parser_m.add_argument("-a","--alpha",default=1,type=float)

parser_flowmap = parser.add_argument_group(description = "==== Flowmap between mixing cycles ====")
parser_flowmap.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")

args=parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

if args.verbose:
    sys.stdout.write(g.ParameterString())
    sys.stdout.write("\n generating matrices\n")

m1,m2 = g.growthmatrixgrid
gm1   = g.growthmatrix[:,:,0]
gm2   = g.growthmatrix[:,:,1]

da = 0.5*(g.growthrates[0] - g.growthrates[1])/(g.growthrates[0] + g.growthrates[1])

da2 = 0.5*da
da2p = 1+da2
da2m = 1-da2


for i,n1 in enumerate(m1):
    for j,n2 in enumerate(m2):
        if n1+n2 > 0:

            step = 0
            n = n1+n2
            x = n1/n
            NN = (gm1[i,j] + gm2[i,j])
            xi0 = NN/n
            
            # set startvalue for iteration
            xi = xi0
            dxi = 1.

            # NR iterations
            while cont(dxi,args.precision,step,args.maxsteps):
                f   = (x * np.power(xi,da2p) + (1-x)*np.power(xi,da2m)) - xi0
                df  = (x * da2p * np.power(xi,da2) + (1-x)*da2m*np.power(xi,-da2))

                dxi = f/df
                xi -= args.alpha * dxi

                step += 1

            print "{:4d} {:4d} {:.6e}".format(n1,n2,xi)
        else:
            print "{:4d} {:4d} {:.6e}".format(n1,n2,0.)
    print ""
                

