#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc


def get_match(n1,ntotal,nlist):
    if ntotal - n1 in nlist:    return ntotal - n1, np.where(nlist == ntotal - n1)[0][0]
    else:                       return None,        None


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",default=None)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-m","--maxM",default=100,type=int)
parser.add_argument("-N","--sortPopSize",default=False,action="store_true")
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError("could not open file '{:s}'".format(args.infile))

openfile = (not args.outfile is None)
if openfile:
    fp = open(args.outfile,"w")
else:
    fp = sys.stdout
    
nlist1,nlist2 = g.growthmatrixgrid
gm1           = g.growthmatrix[:,:,0]
gm2           = g.growthmatrix[:,:,1]
avg_gm1       = np.zeros((len(nlist1),len(nlist2)))
avg_gm2       = np.zeros((len(nlist1),len(nlist2)))

lastn2 = -1

for i,n1 in enumerate(nlist1[nlist1 < args.maxM]):
    for j,n2 in enumerate(nlist2[nlist2 < args.maxM]):
        
        p1 = gc.PoissonSeedingVectors(nlist1,np.array([n1]))
        p2 = gc.PoissonSeedingVectors(nlist2,np.array([n2]))
        
        avg_gm1[i,j] = np.dot(p2[0],np.dot(p1[0],gm1))
        avg_gm2[i,j] = np.dot(p2[0],np.dot(p1[0],gm2))
        
        if not args.sortPopSize:
            
            if lastn2 > n2: fp.write("\n")
            lastn2 = n2
            
            fp.write("{:3.0f} {:3.0f} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n1,n2,avg_gm1[i,j],avg_gm2[i,j],gm1[i,j],gm2[i,j]))

if args.sortPopSize:
    for n in range(min(nlist1[-1] + nlist2[-1],args.maxM)):
        for i,n1 in enumerate(nlist1):

            n2,j = get_match(n1,n,nlist2)

            if not n2 is None:
                fp.write("{:3.0f} {:3.0f} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n1,n2,avg_gm1[i,j],avg_gm2[i,j],gm1[i,j],gm2[i,j]))

        fp.write("\n")

if openfile:
    fp.close()







