#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc
import itertools


parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "Input/Output [required]")
parser_io.add_argument("-i","--infile",required = True)
parser_io.add_argument("-o","--outfile",required = True)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")

parser_dilution = parser.add_argument_group(description = "Parameters for dilution values")
parser_dilution.add_argument("-d","--dilutionmin",type=float,default=1e-6)
parser_dilution.add_argument("-D","--dilutionmax",type=float,default=1e-5)
parser_dilution.add_argument("-K","--dilutionstep",type=float,default=1e-6)

parser_flowmap = parser.add_argument_group(description = "Parameters for Flowmap between mixing cycles")
parser_flowmap.add_argument("-n","--maxIC",type=float,default=40)
parser_flowmap.add_argument("-s","--stepIC",type=float,default=2)
parser_flowmap.add_argument("-L","--trajectorylength",type=int,default=20)

args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

if isinstance(g.growthmatrixgrid,int):
    m = np.arange(g.growthmatrixgrid)
elif isinstance(g.growthmatrixgrid,(tuple,list,np.ndarray)):
    if len(len(g.grothmatrixgrid)) == 1:
        m = g.growthmatrix
    else:
        raise ValueError
else:
    raise ValueError

dlist = np.arange(start = args.dilutionmin,stop = args.dilutionmax + args.dilutionstep,step = args.dilutionstep)
nlist = np.arange(start = 0,stop = args.maxIC,step = args.stepIC)

gm1 = g.growthmatrix[:,:,0]
gm2 = g.growthmatrix[:,:,1]

if args.verbose:
    print g.ParameterString()

for dilution in dlist:
    if args.verbose:
        print "# computing trajectories for D = {:e}".format(dilution)
    fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")
    for icx,icy in itertools.product(nlist,repeat=2):
        x = icx
        y = icy
        print >> fp,x,y
        for i in range(args.trajectorylength):
            px,py = gc.PoissonSeedingVectors(m,np.array((x,y)))
            x = np.dot(py,np.dot(px,gm1))*dilution
            y = np.dot(py,np.dot(px,gm2))*dilution
            print >> fp,x,y
            if (x==0) and (y==0):
                break
        print >> fp
    fp.close()
            
