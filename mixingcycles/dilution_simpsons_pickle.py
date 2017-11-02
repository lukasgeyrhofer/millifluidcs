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


if args.dilutionmax is None:
    dlist = np.array([args.dilutionmin])
else:
    if args.dilutionlogscale:
        dlist = np.power(10,np.linspace(start = np.log10(args.dilutionmin), stop = np.log10(args.dilutionmax), num = args.dilutionbins))
    else:
        dlist = np.linspace(start = args.dilutionmin,stop = args.dilutionmax, num = args.dilutionbins)

nlist = np.arange(start = 0,stop = args.maxIC,step = args.stepIC)
m1,m2 = g.growthmatrixgrid

# growth matrices
gm1 = g.growthmatrix[:,:,0]
gm2 = g.growthmatrix[:,:,1]

# seeding matrices
sm1 = np.array(np.repeat(np.transpose([m1]),len(m2),axis=1),dtype=float)
sm2 = np.repeat(np.array([m2],dtype = float),len(m1),axis=0)

# expressions for simpson's paradoxon
dpp1 = np.nan_to_num(gm1/(gm1 + gm2))
dpp2 = np.nan_to_num(gm2/(gm1 + gm2))

dp1  = np.nan_to_num(sm1/(sm1 + sm2))
dp2  = np.nan_to_num(sm2/(sm1 + sm2))

dps1 = dp2 * gm1
dps2 = dp1 * gm2


if args.verbose:
    print g.ParameterString()

for dilution in dlist:
    if args.verbose:
        print "# computing single step dynamics for D = {:e}".format(dilution)
    fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")
    for x,y in itertools.product(nlist,repeat=2):
        p1 = gc.PoissonSeedingVectors(m1,x)[0]
        p2 = gc.PoissonSeedingVectors(m2,y)[0]
        
        G1 = np.dot(p2,np.dot(p1,gm1)) * dilution
        G2 = np.dot(p2,np.dot(p1,gm2)) * dilution
        
        R1 = np.dot(p2,np.dot(p1,dpp1)) * dilution
        R2 = np.dot(p2,np.dot(p1,dpp2)) * dilution
        
        r1 = np.dot(p2,np.dot(p1,dp1)) * dilution
        r2 = np.dot(p2,np.dot(p1,dp2)) * dilution
        
        S1 = np.dot(p2,np.dot(p1,dps1)) * dilution
        S2 = np.dot(p2,np.dot(p1,dps2)) * dilution
        
        fp.write('{:6.3f} {:6.3f} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n'.format(x,y,G1,G2,R1,R2,r1,r2,S1,S2))
        #                                                                                                           1 2 3  4  5  6  7  8  9  10
            
    fp.close()
            
