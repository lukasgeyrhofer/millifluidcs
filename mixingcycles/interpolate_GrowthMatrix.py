#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc
import itertools

# currently only linear interpolation implemented
def interp_values(w1,w2,edgevalues):
    return (1-w1)*(1-w2)*edgevalues[0] + (1-w1)*w2*edgevalues[1] + w1*(1-w2)*edgevalues[2] + w1*w2*edgevalues[3]

# general helper function that returns strain1/2 inocula, independent of the coordinates used
def mcoord(c1,c2,newcoordinates = False):
    if newcoordinates:
        return np.array([c1*c2,c1*(1-c2)],dtype=np.float)
    else:
        return np.array([c1,c2],dtype=np.float)

parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== Input/Output [required] ====")
parser_io.add_argument("-i","--infile",required = True)
parser_io.add_argument("-o","--outfile",default=None)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")

parser_lattice = parser.add_argument_group(description = "==== New Lattice ====")
parser_lattice.add_argument("-C","--newcoordinates",action="store_true",default=False)
parser_lattice.add_argument("-N","--maxInoculum",type=float,default=100)
parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)
parser_lattice.add_argument("-S","--sorting",choices = ("1","2"),default = "1")
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

mlist1,mlist2 = g.growthmatrixgrid
gm1   = g.growthmatrix[:,:,0]
gm2   = g.growthmatrix[:,:,1]

val1 = np.zeros(4)
val2 = np.zeros(4)

# define new lattice for interpolated growthmatrices
if args.newcoordinates:
    mat1    = gm1 + gm2
    mat2    = np.zeros(np.shape(nmat))
    mat2[mat1 > 0] = gm1[mat1>0]/mat1[mat1>0]

    coordlist1 = np.arange(start = 0,stop = args.maxInoculum, step = args.stepInoculum)
    coordlist2 = np.linspace(start = 0,stop = 1,num = int(1./args.stepFraction)+1)

    new1 = np.zeros((len(coordlist1),len(coordlist2)))
    new2 = np.zeros((len(coordlist1),len(coordlist2)))
else:
    mat1 = gm1
    mat2 = gm2
    
    coordlist1 = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    coordlist2 = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    
    new1  = np.zeros((len(coordlist1),len(coordlist2)))
    new2  = np.zeros((len(coordlist1),len(coordlist2)))

# compute interpolation
for i,c1 in enumerate(coordlist1):
    for j,c2 in enumerate(coordlist2):
        m1i = len(mlist1[mlist1 < mcoord(c1,c2,args.newcoordinates)[0]])
        m2i = len(mlist2[mlist2 < mcoord(c1,c2,args.newcoordinates)[1]])
        if m1i < len(mlist1) and m2i < len(mlist2):
            w1 = (mcoord(c1,c2,args.newcoordinates)[0] - mlist1[m1i])/float(mlist1[m1i+1] - mlist1[m1i])
            w2 = (mcoord(c1,c2,args.newcoordinates)[1] - mlist2[m2i])/float(mlist2[m2i+1] - mlist2[m2i])

            val1 = np.array([nmat[m1i, m2i], nmat[m1i, m2i+1], nmat[m1i+1, m2i], nmat[m1i+1, m2i+1]])
            val2 = np.array([xmat[m1i, m2i], xmat[m1i, m2i+1], xmat[m1i+1, m2i], xmat[m1i+1, m2i+1]])

            new1[i,j] = interp_values(w1,w2,valn)
            new2[i,j] = interp_values(w1,w2,valx)


# output
if not args.outfile is None:
    fp = open(args.outfile,"w")
else:
    fp = sys.stdout

if args.sorting == '1':
    for i,c1 in enumerate(coordlist1):
        for j,c2 in enumerate(coordlist2):
            fp.write("{:7.2f} {:.6f} {:.6e} {:.6e}\n".format(c1,c2,new1[i,j],new2[i,j]))
        fp.write("\n")
else:
    for j,c2 in enumerate(coordlist2):
        for i,c1 in enumerate(coordlist1):
            fp.write("{:7.2f} {:.6f} {:.6e} {:.6e}\n".format(c1,c2,new1[i,j],new2[i,j]))
        fp.write("\n")

if not args.outfile is None:
    fp.close()



