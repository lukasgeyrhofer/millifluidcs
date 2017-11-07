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

print "=================================================================="
print "  growthmatrix file : {}".format(args.infile)
print "  dynamics type :     {}".format(str(type(g)).split("'")[1])
print "  GMgrid :            {}".format(g.growthmatrixgrid)
print "  GMshape :           {}".format(np.shape(g.growthmatrix))
print "=================================================================="
print

print g.ParameterString()
