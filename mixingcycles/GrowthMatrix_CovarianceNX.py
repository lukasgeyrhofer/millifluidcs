#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc

def main():
    parser = argparse.ArgumentParser()

    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infile")
    parser_io.add_argument("-o","--outfile",default=None)

    parser_lattice = parser.add_argument_group(description = "==== Lattice parameters ====")
    parser_lattice.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
    parser_lattice.add_argument("-N","--maxInoculum",type=float,default=40)
    parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
    parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)

    args = parser.parse_args()

    
    # load growthmatrix and growthmatrix grid
    try:
        g = pickle.load(open(args.infile))
    except:
        raise IOError,"Could not open and load from pickle file"

    if not g.hasGrowthMatrix():
        raise ValueError,"Loaded pickle instance does not contain growthmatrix"

    gm1          = g.growthmatrix[:,:,0]
    gm2          = g.growthmatrix[:,:,1]
    m1,m2        = g.growthmatrixgrid

    # transformations
    # initial conditions as matrices
    mm1          = np.repeat([m1.T],len(m2),axis=0).T
    mm2          = np.repeat([m2],  len(m1),axis=0)
    nIC          = 1.*mm1 + mm2
    xIC          = np.zeros(np.shape(mm1))
    xIC[nIC > 0] = 1.*mm1[nIC>0]/nIC[nIC>0]

    # final sizes as matrices
    nn           = gm1 + gm2
    xx           = np.zeros(np.shape(nn))
    xx[nn>0]     = gm1[nn>0]/nn[nn>0]


    # define new grid
    if args.newcoordinates:
        nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
        xlist = np.arange(start = 0,stop = 1 + .5*args.stepFraction,step = args.stepFraction)
        shape = (len(nlist),len(xlist))
    else:
        nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
        xlist = None
        shape = (len(nlist),len(nlist))

    # observables on new grid
    ENX          = np.zeros(shape)
    EX           = np.zeros(shape)
    EN           = np.zeros(shape)
    Cov          = np.zeros(shape)
    CovNorm      = np.zeros(shape)
    EdX          = np.zeros(shape)

    # output definitions
    if args.outfile is None:
        fpout = sys.stdout
    else:
        fpout = open(args.outfile,"w")


    # iterate over all possible new coodinates
    if args.newcoordinates:
        for i,n in enumerate(nlist):
            for j,x in enumerate(xlist):
                p1 = gc.PoissonSeedingVectors(m1,[n*x])[0]
                p2 = gc.PoissonSeedingVectors(m2,[n*(1-x)])[0]
                
                EdX[i,j] = np.dot(p2, np.dot(xx - xIC, p1))
                ENX[i,j] = np.dot(p2, np.dot(nn * xx,  p1))
                EN[i,j]  = np.dot(p2, np.dot(nn     ,  p1))
                EX[i,j]  = np.dot(p2, np.dot(xx     ,  p1))
                
                Cov[i,j] = (ENX[i,j] - EN[i,j] * EX[i,j])
                CovNorm[i,j] = Cov[i,j]
                relN = 0
                if EN[i,j] > 0:
                    CovNorm[i,j] /= EN[i,j]
                    relN = n/EN[i,j]
               
                fpout.write("{:.6e} {:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n,x,CovNorm[i,j],Cov[i,j],ENX[i,j],EX[i,j],EN[i,j],EdX[i,j],relN))
            fpout.write("\n")
    else:
        for i,n1 in enumerate(nlist):
            for j,n2 in enumerate(nlist):
                p1 = gc.PoissonSeedingVectors(m1,[n1])[0]
                p2 = gc.PoissonSeedingVectors(m2,[n2])[0]
                
                x = 0
                if n1 + n2 > 0:
                    x = 1.*n1/(1.*n1 + n2)
                
                EdX[i,j] = np.dot(p2, np.dot(xx - xIC, p1))
                ENX[i,j] = np.dot(p2, np.dot(nn * xx,  p1))
                EN[i,j]  = np.dot(p2, np.dot(nn     ,  p1))
                EX[i,j]  = np.dot(p2, np.dot(xx     ,  p1))
                
                Cov[i,j] = (ENX[i,j] - EN[i,j] * EX[i,j])
                CovNorm[i,j] = Cov[i,j]
                relN = 0
                if EN[i,j] > 0:
                    CovNorm[i,j] /= EN[i,j]
                    relN = (n1+n2)/EN[i,j]
                
                fpout.write("{:.6e} {:.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n1,n2,CovNorm[i,j],Cov[i,j],ENX[i,j],EX[i,j],EN[i,j],EdX[i,j],relN))
            fpout.write("\n")


if __name__ == "__main__":
    main()
    

