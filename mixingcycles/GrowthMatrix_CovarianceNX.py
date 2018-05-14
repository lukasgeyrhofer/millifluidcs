#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc

parser = argparse.ArgumentParser()

parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
parser_io.add_argument("-i","--infile")
parser_io.add_argument("-o","--outfile",default=None)

parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
parser_alg.add_argument("-L","--normalize",default=False,action="store_true",help="compute Cov[N,X]/E[N] instead of Cov[N,X]")

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

mx,my = g.growthmatrixgrid

if args.newcoordinates:
    nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    xlist = np.arange(start = 0,stop = 1 + .5*args.stepFraction,step = args.stepFraction)
    shape = (len(nlist),len(xlist))
else:
    nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    xlist = None
    shape = (len(nlist),len(nlist))


m1,m2    = g.growthmatrixgrid
gm1      = g.growthmatrix[:,:,0]
gm2      = g.growthmatrix[:,:,1]

nn       = gm1 + gm2
xx       = np.zeros(np.shape(nn))
xx[nn>0] = gm1[nn>0]/nn[nn>0]

Enx      = np.zeros(shape)
Ex       = np.zeros(shape)
En       = np.zeros(shape)
cov      = np.zeros(shape)
covN     = np.zeros(shape)

if args.outfile is None:
    fpout = sys.stdout
else:
    fpout = open(args.outfile,"w")

if args.newcoordinates:
    for i,n in enumerate(nlist):
        for j,x in enumerate(xlist):
            p1 = gc.PoissonSeedingVectors(m1,[n*x])[0]
            p2 = gc.PoissonSeedingVectors(m1,[n*(1-x)])[0]
            
            Enx[i,j] = np.dot(p2, np.dot(nn * xx, p1))
            En[i,j]  = np.dot(p2, np.dot(nn     , p1))
            Ex[i,j]  = np.dot(p2, np.dot(xx     , p1))
            
            cov[i,j] = (Enx[i,j] - En[i,j] * Ex[i,j])
            covN[i,j] = cov[i,j]
            if En[i,j] > 0: covN[i,j] /= En[i,j]
            
            fpout.write("{:.6e} {:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n,x,cov[i,j],covN[i,j],Enx[i,j],Ex[i,j],En[i,j]))
        fpout.write("\n")
else:
    for i,n1 in enumerate(nlist):
        for j,n2 in enumerate(nlist):
            p1 = gc.PoissonSeedingVectors(m1,[n1])[0]
            p2 = gc.PoissonSeedingVectors(m1,[n2])[0]
            
            Enx[i,j] = np.dot(p2, np.dot(nn * xx, p1))
            En[i,j]  = np.dot(p2, np.dot(nn     , p1))
            Ex[i,j]  = np.dot(p2, np.dot(xx     , p1))

            cov[i,j] = (Enx[i,j] - En[i,j] * Ex[i,j])
            covN[i,j] = cov[i,j]
            if En[i,j] > 0: covN[i,j] /= En[i,j]
            
            fpout.write("{:.6e} {:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n,x,cov[i,j],covN[i,j],Enx[i,j],Ex[i,j],En[i,j]))
        fpout.write("\n")


        
