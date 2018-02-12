#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc
import itertools


parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== Input/Output [required] ====")
parser_io.add_argument("-i","--infile",required = True)
parser_io.add_argument("-o","--outfile",required = True)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")

parser_dilution = parser.add_argument_group(description = "==== Parameters for dilution values ====")
parser_dilution.add_argument("-d","--dilutionmin",type=float,default=1e-6)
parser_dilution.add_argument("-D","--dilutionmax",type=float,default=None)
parser_dilution.add_argument("-K","--dilutionsteps",type=int,default=10)
parser_dilution.add_argument("-L","--dilutionlogscale",default = False, action = "store_true")

parser_flowmap = parser.add_argument_group(description = "==== Flowmap between mixing cycles ====")
parser_flowmap.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
parser_flowmap.add_argument("-N","--maxInoculum",type=float,default=40)
parser_flowmap.add_argument("-n","--stepInoculum",type=float,default=2)
parser_flowmap.add_argument("-x","--stepFraction",type=float,default=.05)

args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

mx,my = g.growthmatrixgrid
if args.dilutionmax is None:
    dlist = np.array([args.dilutionmin])
else:
    if args.dilutionlogscale:
        dlist = np.power(10,np.linspace(start = np.log10(args.dilutionmin),stop = np.log10(args.dilutionmax), num = args.dilutionsteps))
    else:
        dlist = np.linspace(start = args.dilutionmin,stop = args.dilutionmax,num = args.dilutionsteps)

if args.newcoordinates:
    nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    xlist = np.arange(start = 0,stop = 1 + .5*args.stepFraction,step = args.stepFraction)
    shape = (len(nlist),len(xlist))
else:
    nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    xlist = None
    shape = (len(nlist),len(nlist))

gm1 = g.growthmatrix[:,:,0]
gm2 = g.growthmatrix[:,:,1]

if args.verbose:
    sys.stdout.write(g.ParameterString())

for dilution in dlist:
    if args.verbose:
        sys.stdout.write("# computing single step dynamics for D = {:e}\n".format(dilution))
    fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")

    if not args.newcoordinates:
        lastx = None
        for x,y in itertools.product(nlist,repeat=2):
            px = gc.PoissonSeedingVectors(mx,x)[0]
            py = gc.PoissonSeedingVectors(my,y)[0]
            nx = np.dot(py,np.dot(px,gm1))*dilution
            ny = np.dot(py,np.dot(px,gm2))*dilution
            if x != lastx:
                fp.write("\n")
            lastx = x
            fp.write('{} {} {} {}\n'.format(x,y,nx,ny))
    else:
        lastx = None
        for x,n in itertools.product(xlist,nlist):
            p1 = gc.PoissonSeedingVectors(mx,n*x)[0]
            p2 = gc.PoissonSeedingVectors(my,n*(1-x))[0]
            n1 = np.dot(p2,np.dot(p1,gm1))*dilution
            n2 = np.dot(p2,np.dot(p1,gm2))*dilution
            nn = n1 + n2
            if nn > 0:  xx = n1/nn
            else:       xx = 0
            if x != lastx:
                fp.write("\n")
            fp.write('{} {} {} {}\n'.format(n,x,nn,xx))

    fp.close()
            
