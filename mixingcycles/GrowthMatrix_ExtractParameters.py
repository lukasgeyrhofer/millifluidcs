#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle
sys.path.append(sys.path[0] + '/../')
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
print "  GMshape :           {}".format(np.shape(g.growthmatrix))
print "  GMgridX :           \n{}".format(g.growthmatrixgrid[0])
print "  GMgridY :           \n{}".format(g.growthmatrixgrid[1])
print "=================================================================="
print

print g.ParameterString()
