#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/..')

import growthclasses as gc
from skimage import measure

parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-i","--infile",required=True)
parser_io.add_argument("-o","--baseoutfilename",default="out")
parser_io.add_argument("-v","--verbose",action="store_true",default=False)


parser_lattice = parser.add_argument_group(description = "==== Lattice parameters ====")
parser_lattice.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
parser_lattice.add_argument("-N","--maxInoculum",type=float,default=40)
parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)

parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
parser_alg.add_argument("-T","--thresholds",nargs="*",type=float)

args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError

m1,m2 = g.growthmatrixgrid
gm1   = g.growthmatrix[:,:,0]
gm2   = g.growthmatrix[:,:,1]

axis1,axis2 = gc.getInoculumAxes(**vars(args))
avgPopSize  = np.zeros((len(axis1),len(axis2)))

for i,a in enumerate(axis1):
    for j,b in enumerate(axis2):
        avgPopSize[i,j] = gc.SeedingAverage(gm1 + gm2,gc.getAbsoluteInoculumNumbers([a,b],args.newcoordinates))
        
avgPopSize /= np.max(avgPopSize)

if len(args.thresholds) > 0:
    for t in args.thresholds:
        fp = open(args.baseoutfilename + '_{:.6f}'.format(t),'w')
        contours = measure.find_contours(avgPopSize,t)
        for c in contours:
            out2 = list()
            for i in range(len(c)):
                ix = int(np.floor(c[i,0]))
                iy = int(np.floor(c[i,1]))
                px = c[i,0] - ix
                py = c[i,1] - iy
                
                try:    cx = (1.-px)*axis1[ix] + px*axis1[ix+1]
                except: cx = axis1[ix]
                try:    cy = (1.-py)*axis2[iy] + py*axis2[iy+1]
                except: cy = axis2[iy]
        
                fp.write('{} {}\n'.format(cx,cy))
            fp.write('\n')
        fp.close()
    
    

