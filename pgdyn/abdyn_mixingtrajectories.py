#!/usr/bin/env python

import argparse
import numpy as np
import itertools
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=False)

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
parser_iterationmap.add_argument("-N","--outputdx",type=float,default=.1)
#parser_iterationmap.add_argument("-P","--poissonseeding",default=False,action="store_true")
parser_iterationmap.add_argument("-t","--trajectorylength",type=int,default=20)

parser_dilution = parser.add_argument_group(description = "Parameters for dilution values")
parser_dilution.add_argument("-d","--dilutionmin",type=float,default=1e-6)
parser_dilution.add_argument("-D","--dilutionmax",type=float,default=1e-5)
parser_dilution.add_argument("-K","--dilutionstep",type=float,default=1e-6)

parser.add_argument("-o","--outfile",default=None)


args      = parser.parse_args()
g         = gc.GrowthDynamicsAntibiotics(**vars(args))
gm1,gm2   = g.getGrowthMatrix(size = args.maxsize,step = args.step)
m         = np.arange(start = 0,stop = args.maxsize,step = args.step,dtype=int)
outpoints = np.arange(start = 0,stop = args.outputmax,step = args.outputdx,dtype=float)
dilutions = np.arange(start = args.dilutionmin,stop = args.dilutionmax,step = args.dilutionstep,dtype=float)


if args.outfile is None:
    print >> sys.stderr,"need outfile"
    exit(1)

for dilution in dilutions:
    outfile = "%s_D%f"%(args.outfile.replace(" ",""),dilution)
    fp = open(outfile,"w")
    
    for icx,icy in itertools.product(outpoints,repeat=2)
        x = icx
        y = icy
        print >> fp,x,y
        for i in range(args.trajectorylength):
            px,py = gc.PoissonSeedingVectors(m,np.array((x,y)))
            x = np.dot(py,np.dot(px,gm1))
            y = np.dot(py,np.dot(px,gm2))
            print >> fp,x,y
        print >> fp
    
    fp.close()
