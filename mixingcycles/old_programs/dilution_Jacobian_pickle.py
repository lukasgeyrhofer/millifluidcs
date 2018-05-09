#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-d","--dilution",type=float,default=1e-4)
parser.add_argument("-n","--inoculumsize",nargs="*",type=float,default=[1.,1.])
args = parser.parse_args()

n = np.array(args.inoculumsize,dtype=np.float)

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"could not open pickle file"


if not g.hasGrowthMatrix():
    raise ValueError,"pickle file does not contain growthmatrix"

gm1 = g.growthmatrix[:,:,0]
gm2 = g.growthmatrix[:,:,1]

m = np.arange(g.growthmatrixgrid)

px,dpx = gc.PoissonSeedingVectors(m,n,diff=True)

jac = np.array([[ np.dot( px[1],np.dot(dpx[0],gm1))*args.dilution ,
                  np.dot(dpx[1],np.dot( px[0],gm1))*args.dilution ],
                [ np.dot( px[1],np.dot(dpx[0],gm2))*args.dilution ,
                  np.dot(dpx[1],np.dot( px[0],gm2))*args.dilution ]],dtype=float)

eigval,eigvec = np.linalg.eig(jac)

print jac

print eigval
print eigvec[0]
print eigvec[1]

