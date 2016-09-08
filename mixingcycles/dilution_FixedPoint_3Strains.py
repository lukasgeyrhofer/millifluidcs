#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultgrowthrates = [1.,1.,1.],defaultyieldfactors = [1.,1.,1.],defaultdeathrates = [0.,0.,0.],dilution = True)
parser.add_argument("-m","--maxseeding",type=int,default=200)
parser.add_argument("-p","--precision",type=float,default=1e-10)
parser.add_argument("-A","--alpha",type=float,default=1)
parser.add_argument("-M","--maxiterations",type=int,default=10000)
parser.add_argument("-N","--newtonraphson",action="store_true",default=False)
parser.add_argument("-v","--verbose",action="store_true",default=False)
args = parser.parse_args()

g = gc.GrowthDynamics(**vars(args))

assert g.numstrains == 3

print("# creating growth matrices ...")
g1,g2,g3 = g.getGrowthMultipleStrains(size = args.maxseeding,nstrains = g.numstrains)
print("# growth matrices done ...")

m = np.arange(args.maxseeding)
J = np.zeros((3,3))
i=0

n  = np.ones(g.numstrains) * g.env.dilution * np.mean(g.yieldfactors) * g.env.substrate
dn = n

while np.sum((dn[n>0]/n[n>0])**2) > args.precision:
    px,dpx = gc.PoissonSeedingVectors(m,n,diff=True)
    
    growth1 = np.dot(px[2],np.dot(px[1],np.dot(px[0],g1)))
    growth2 = np.dot(px[2],np.dot(px[1],np.dot(px[0],g2)))
    growth3 = np.dot(px[2],np.dot(px[1],np.dot(px[0],g3)))
    
    fn = np.array([growth1,growth2,growth2]) - n
    
    if args.newtonraphson:
        
        J[0,0] = np.dot( px[2], np.dot( px[1], np.dot(dpx[0],g1)))-1.
        J[0,1] = np.dot( px[2], np.dot(dpx[1], np.dot( px[0],g1)))
        J[0,2] = np.dot(dpx[2], np.dot( px[1], np.dot( px[0],g1)))
                                                       
        J[1,0] = np.dot( px[2], np.dot( px[1], np.dot(dpx[0],g1)))
        J[1,1] = np.dot( px[2], np.dot(dpx[1], np.dot( px[0],g1)))-1.
        J[1,2] = np.dot(dpx[2], np.dot( px[1], np.dot( px[0],g1)))
                                                       
        J[2,0] = np.dot( px[2], np.dot( px[1], np.dot(dpx[0],g1)))
        J[2,1] = np.dot( px[2], np.dot(dpx[1], np.dot( px[0],g1)))
        J[2,2] = np.dot(dpx[2], np.dot( px[1], np.dot( px[0],g1)))-1.
        
        dn = -args.alpha * np.dot(np.linalg.inv(J),fn)
    
    else:
        dn = fn
    
    n += dn
    n[n<0] = 0.

    if not args.maxiterations is None:
        if i > args.maxiterations:
            break
    i += 1


    print("{:5d} {:10.6f} {:10.6f} {:10.6f}".format(i,n[0],n[1],n[2]))
