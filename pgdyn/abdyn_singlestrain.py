#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


sys.path.append(sys.path[0] + '/../mixingcyclyes')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parserAB = parser.add_argument_group(description = "Parameters for interactions with antibiotics")
parserAB.add_argument("-k","--kappa",type=float,default=1)
parserAB.add_argument("-l","--logkill",type=float,default=2)
parserAB.add_argument("-P","--PGproduction",nargs="*")
parserAB.add_argument("-A","--ABconc",type=float,default=.5)
parserAB.add_argument("-R","--PGreductionAB",type=float,default=1)

parser.add_argument("-m","--maxsize",type=int,default=100)
parser.add_argument("-n","--outputmax",type=float,default=10)
parser.add_argument("-d","--outputdx",type=float,default=.1)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-p","--poissonseeding",default=False,action="store_true")

args = parser.parse_args()

g = gc.GrowthDynamicsAntibiotics(**vars(args))
m = np.arange(args.maxsize)

gm1 = g.getGrowthVector(args.maxsize)

if args.poissonseeding:
    out    = np.arange(0,args.outputmax+args.outputdx,args.outputdx)
    p      = gc.PoissonSeedingVectors(m,out)
    growth = np.dot(p,gm1)
else:
    out = m[:]
    growth = gm1[:]
for a,b in zip(out,growth):
    print a,b

