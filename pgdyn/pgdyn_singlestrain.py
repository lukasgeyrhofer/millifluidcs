#!/usr/bin/env python

import numpy as np
import argparse
import sys

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser_PG = parser.add_argument_group(description = "Parameters for interactions with public good")
parser_PG.add_argument("-P","--pgproduction",nargs="*",default=np.zeros(2))
parser_PG.add_argument("-A","--pginteractiongrowthrates",nargs="*",default=np.zeros(2))
parser_PG.add_argument("-Y","--pginteractionyieldfactor",nargs="*",default=np.zeros(2))
parser_PG.add_argument("-O","--onlypositivecoefficiets",action="store_true",default=False)

parser.add_argument("-m","--maxsize",type=int,default=100)
parser.add_argument("-n","--outputmax",type=float,default=10)
parser.add_argument("-d","--outputdx",type=float,default=.1)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-p","--poissonseeding",default=False,action="store_true")

args = parser.parse_args()

g = gc.GrowthDynamicsPublicGoods(**vars(args))
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




