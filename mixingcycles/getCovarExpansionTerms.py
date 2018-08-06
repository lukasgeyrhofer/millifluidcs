#!/usr/bin/env python


import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc


parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-i","--infile",required=True)
parser_io.add_argument("-o","--baseoutfilename",default="out")
parser_io.add_argument("-v","--verbose",action="store_true",default=False)


parser_lattice = parser.add_argument_group(description = "==== Lattice parameters ====")
parser_lattice.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
parser_lattice.add_argument("-N","--maxInoculum",type=float,default=40)
parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)

args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

if args.verbose:
    sys.stdout.write(g.ParameterString())

if args.newcoordinates:
    nlist = np.arange(start = 0,stop = args.maxInoculum + .5*args.stepInoculum,step = args.stepInoculum)
    xlist = np.arange(start = 0,stop = 1 + .5*args.stepFraction,step = args.stepFraction)
    shape = (len(nlist),len(xlist))
else:
    nlist = np.arange(start = 0,stop = args.maxInoculum+ .5*args.stepInoculum,step = args.stepInoculum)
    xlist = None
    shape = (len(nlist),len(nlist))


m1,m2 = g.growthmatrixgrid
gm1   = g.growthmatrix[:,:,0]
gm2   = g.growthmatrix[:,:,1]

mm1   = np.repeat([m1],len(m2),axis=0)
mm2   = np.repeat([m2],len(m1),axis=0).T

wd_Nfin   = np.zeros(gm1.shape)
wd_Nini   = np.zeros(gm1.shape)
wd_Xfin   = np.zeros(gm1.shape)
wd_Xini   = np.zeros(gm1.shape)

wd_Nfin            = gm1 + gm2
wd_Nini            = mm1 + mm2
wd_Xfin[wd_Nfin>0] = gm1[wd_Nfin>0]/wd_Nfin[wd_Nfin>0]
wd_Xini[wd_Nini>0] = mm1[wd_Nini>0]/wd_Nini[wd_Nini>0]
wd_dX              = wd_Xfin - wd_Xini
wd_XNNfin          = gm1 * wd_Nfin

wd_xi = g.GetXiMatrix()

if args.newcoordinates:
    for i,n in enumerate(nlist):
        for j,x in enumerate(xlist):
            
            avg_Nfin  = gc.SeedingAverage(wd_Nfin,         n*x, n*(1-x))
            avg_N1fin = gc.SeedingAverage(gm1,             n*x, n*(1-x))
            avg_Xfinwd= gc.SeedingAverage(wd_Xfin,         n*x, n*(1-x))
            avg_N1N   = gc.SeedingAverage(gm1 * wd_Nfin,   n*x, n*(1-x))
            avg_nxi   = gc.SeedingAverage(wd_Nini * wd_xi, n*x, n*(1-x))
            
            avg_Xfin  = avg_N1fin / avg_Nfin
            
            cov       = (avg_N1N - avg_N1fin * avg_Nfin)/(avg_Nfin * avg_Nfin)
            Xomega1   = wd_Xfin * (wd_Nini * wd_xi / avg_nxi - 1)
            
            avg_Xomega1 = gc.SeedingAverage(Xomega1,n*x, n*(1-x))

            print "{:10.2f} {:.6f} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(n,x,avg_Nfin,cov,avg_Xomega1,avg_Xfin,avg_Xfinwd)
        print
else:
    for i,n1 in enumerate(nlist):
        for j,n2 in enumerate(nlist):
            
            avg_Nfin  = gc.SeedingAverage(wd_Nfin,         n1, n2)
            avg_N1fin = gc.SeedingAverage(gm1,             n1, n2)
            avg_Xfinwd= gc.SeedingAverage(wd_Xfin,         n1, n2)
            avg_N1N   = gc.SeedingAverage(gm1 * wd_Nfin,   n1, n2)
            avg_nxi   = gc.SeedingAverage(wd_Nini * wd_xi, n1, n2)
            
            avg_Xfin  = avg_N1fin / avg_Nfin
            
            cov       = (avg_N1N - avg_N1fin * avg_Nfin)/(avg_Nfin * avg_Nfin)
            Xomega1   = X * (wd_Nini * wd_xi / avg_nxi - 1)
            
            avg_Xomega1 = gc.SeedingAverage(Xomega1,n1, n2)

            print "{:10.2f} {:.6f} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(n,x,avg_Nfin,cov,avg_Xomega1,avg_Xfin,avg_Xfinwd)
        print
            
            
            








