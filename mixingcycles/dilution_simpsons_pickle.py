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
parser_dilution.add_argument("-D","--dilutionmax",type=float,default=None)
parser_dilution.add_argument("-K","--dilutionbins",type=int,default=20)
parser_dilution.add_argument("-L","--dilutionlogscale",default=False,action="store_true")

parser_flowmap = parser.add_argument_group(description = "Parameters for Flowmap between mixing cycles")
parser_flowmap.add_argument("-n","--maxIC",type=float,default=40)
parser_flowmap.add_argument("-s","--stepIC",type=float,default=2)

args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

m1,m2 = g.growthmatrixgrid

m1A = np.repeat([m1],len(m2),axis=1)
m2A = np.repeat(np.transpose([m2]),len(m1),axis=0)

print m1A
print m1A[0,2]
exit(1)

if args.dilutionmax is None:
    dlist = np.array([args.dilutionmin])
else:
    if args.dilutionlogscale:
        dlist = np.power(10,np.linspace(start = np.log10(args.dilutionmin), stop = np.log10(args.dilutionmax), num = args.dilutionbins))
    else:
        dlist = np.linspace(start = args.dilutionmin,stop = args.dilutionmax, num = args.dilutionbins)

nlist = np.arange(start = 0,stop = args.maxIC,step = args.stepIC)

gm1 = g.growthmatrix[:,:,0]
gm2 = g.growthmatrix[:,:,1]

if args.verbose:
    print g.ParameterString()

for dilution in dlist:
    if args.verbose:
        print "# computing single step dynamics for D = {:e}".format(dilution)
    fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")
    for x,y in itertools.product(nlist,repeat=2):
        px = gc.PoissonSeedingVectors(mx,x)[0]
        py = gc.PoissonSeedingVectors(my,y)[0]
        nx = np.dot(py,np.dot(px,gm1))*dilution
        ny = np.dot(py,np.dot(px,gm2))*dilution
        fp.write('{} {} {} {}'.format(x,y,nx,ny))
            
    fp.close()
            
