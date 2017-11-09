#!/usr/bin/env python

import argparse
import numpy as np
import sys
sys.path.append(sys.path[0] + '/../')
import growthclasses as gc
import pickle
import inspect


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile1")
parser.add_argument("-I","--infile2")
parser.add_argument("-o","--outfile")
args = parser.parse_args()

try:
    g1 = pickle.load(open(args.infile1))
    g2 = pickle.load(open(args.infile2))
except:
    raise IOError("could not open files")

gm = g1._GrowthDynamics__growthmatrix
g2._GrowthDynamics__growthmatrix[:len(gm[:,0,0]),:len(gm[0,:,0]),:] = gm

fp = open(args.outfile,"w")
pickle.dump(g2,fp)
fp.close()

