#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",default=None)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-m","--maxM",default=100)
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError("could not open file '{:s}'".format(args.infile))

openfile = (not args.outfile is None)
if openfile:
    fp = open(args.outfile)
else:
    fp = sys.stdout
    
nlist1  = g.__growthmatrixgridX[g.__growthmatrixgridX < args.maxM]
nlist2  = g.__growthmatrixgridY[g.__growthmatrixgridY < args.maxM]
gm1     = g.__growthmatrix[::,0]
gm2     = g.__growthmatrix[::,1]
avg_gm1 = np.zeros((len(nlist1),len(nlist2)))
avg_gm2 = np.zeros((len(nlist1),len(nlist2)))

lastn2 = -1
for i,n1 in enumerate(nlist1):
    for j,n2 in enumerate(nlist2):
        
        p1 = gc.PoissonSeedingVectors(g.__growthmatrixgridX,[n1])
        p2 = gc.PoissonSeedingVectors(g.__growthmatrixgridY,[n2])
        
        avg_gm1[i,j] = np.dot(p2[0],np.dot(p1[0],gm1))
        avg_gm2[i,j] = np.dot(p2[0],np.dot(p1[0],gm2))
        
        if lastn2 > n2: fp.write("\n")
        lastn2 = n2
        
        fp.write("{:3.0f} {:3.0f} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(n1,n2,avg_gm1[i,j],avg_gm2[i,j],gm1[i,j],gm2[i,j]))

if openfile:
    fp.close()







