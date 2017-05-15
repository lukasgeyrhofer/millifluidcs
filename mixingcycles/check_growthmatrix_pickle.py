#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",default=None)
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

print g.ParameterString()

print "*** GrowthMatrixGrid ***"
print "  GMgrid   {}".format(g.growthmatrixgrid)
print "  GMshape  {}".format(np.shape(g.growthmatrix))
