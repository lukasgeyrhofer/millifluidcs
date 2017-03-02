#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=True)

parser_ab = parser.add_argument_group(description = "Parameters for dynamics of AB")
parser_ab.add_argument("-B","--ABconc",type=float,default=1.2)
parser_ab.add_argument("-g","--gamma",type=float,default=2)
parser_ab.add_argument("-k","--kappa",type=float,default=2)
parser_ab.add_argument("-p","--PGproduction",nargs="*",default=[1,0])
parser_ab.add_argument("-r","--PGreductionAB",type=float,default=1e-3)

parser_iterationmap = parser.add_argument_group(description = "Parameters for iterationmap")
parser_iterationmap.add_argument("-m","--maxsize",type=int,default=100)
parser_iterationmap.add_argument("-M","--step",type=int,default=1)
parser_iterationmap.add_argument("-n","--outputmax",type=float,default=20)
parser_iterationmap.add_argument("-d","--outputdx",type=float,default=.1)
parser_iterationmap.add_argument("-P","--poissonseeding",default=False,action="store_true")

parser.add_argument("-o","--outfile",default=None)

args = parser.parse_args()
g = gc.GrowthDynamicsAntibiotics(**vars(args))
gm1,gm2 = g.getGrowthMatrix(size = args.maxsize)

m = np.arange(args.maxsize)*args.step
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
        print >> fp,"{} {} {} {}".format(x,y,gx,gy)
    print >> fp

if not args.outfile is None:
    fp.close()
