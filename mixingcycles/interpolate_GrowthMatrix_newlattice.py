#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc
import itertools


def interp_values(w1,w2,edgevalues):
    return (1-w1)*(1-w2)*edgevalues[0] + (1-w1)*w2*edgevalues[1] + w1*(1-w2)*edgevalues[2] + w1*w2*edgevalues[3]

parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== Input/Output [required] ====")
parser_io.add_argument("-i","--infile",required = True)
parser_io.add_argument("-o","--outfile",default=None)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")

parser_lattice = parser.add_argument_group(description = "==== Lattice in transformed coordinates ====")
parser_lattice.add_argument("-N","--maxInoculum",type=float,default=40)
parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)
parser_lattice.add_argument("-S","--sorting",choices = ("n","x"),default = "n")
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

mx,my = g.growthmatrixgrid
gm1   = g.growthmatrix[:,:,0]
gm2   = g.growthmatrix[:,:,1]

nmat    = gm1 + gm2
xmat    = np.zeros(np.shape(nmat))
xmat[nmat > 0] = gm1[nmat>0]/nmat[nmat>0]

nlist = np.arange(start = 0,stop = args.maxInoculum, step = args.stepInoculum)
xlist = np.linspace(start = 0,stop = 1,num = int(1./args.stepFraction)+1)

newn = np.zeros((len(nlist),len(xlist)))
newx = np.zeros((len(nlist),len(xlist)))
valn = np.zeros(4)
valx = np.zeros(4)


for i,n in enumerate(nlist):
    for j,x in enumerate(xlist):
        m1 = int(np.floor(n*x))
        m2 = int(np.floor(n*(1-x)))
        w1 = n*x - m1
        w2 = n*(1-x) - m2
        if m1 < mx[-1] and m2 < my[-1]:
            edges = np.array([[m1,m2],[m1,m2+1],[m1+1,m2],[m1+1,m2+1]])
            valn[0] = nmat[m1,   m2]
            valn[1] = nmat[m1,   m2+1]
            valn[2] = nmat[m1+1, m2]
            valn[3] = nmat[m1+1, m2+1]

            valx[0] = xmat[m1,   m2]
            valx[1] = xmat[m1,   m2+1]
            valx[2] = xmat[m1+1, m2]
            valx[3] = xmat[m1+1, m2+1]

            newn[i,j] = interp_values(w1,w2,valn)
            newx[i,j] = interp_values(w1,w2,valx)

if not args.outfile is None:
    fp = open(args.outfile,"w")
else:
    fp = sys.stdout

if args.sorting == 'n':
    for i,n in enumerate(nlist):
        for j,x in enumerate(xlist):
            fp.write("{:7.2f} {:.6f} {:.6e} {:.6e}\n".format(n,x,newn[i,j],newx[i,j]))
        fp.write("\n")
else:
    for j,x in enumerate(xlist):
        for i,n in enumerate(nlist):
            fp.write("{:7.2f} {:.6f} {:.6e} {:.6e}\n".format(n,x,newn[i,j],newx[i,j]))
        fp.write("\n")

if not args.outfile is None:
    fp.close()



