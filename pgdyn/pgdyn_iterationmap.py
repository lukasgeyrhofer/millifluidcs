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

parser.add_argument("-m","--maxsize",type=int,default=100)
parser.add_argument("-n","--outputmax",type=float,default=10)
parser.add_argument("-d","--outputdx",type=float,default=.1)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-P","--poissonseeding",default=False,action="store_true")

args = parser.parse_args()

g       = gc.GrowthDynamicsPublicGoods(**vars(args))
gm1,gm2 = g.getGrowthMatrix(size = args.maxsize)

m = np.arange(args.maxsize)
if args.poissonseeding:
    outpoints = np.arange(0,args.outputmax,args.outputdx,dtype=float)
else:
    outpoints = np.arange(args.maxsize,dtype=int)
    

writetofile = False
if not args.outfile is None:
    try:
        fp = open(args.outfile,"w")
        writetofile = True
    except:
        fp = sys.stdout
else:
    fp = sys.stdout
    

for x in outpoints:
    for y in outpoints:
        if args.poissonseeding:
            px,py = gc.PoissonSeedingVectors(m,np.array((x,y)))
            gx = np.dot(py,np.dot(px,gm1))
            gy = np.dot(py,np.dot(px,gm2))
        else:
            gx = gm1[x,y]
            gy = gm2[x,y]
        print "{} {} {} {}".format(x,y,gx,gy)
    print

