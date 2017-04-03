#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile",default=None)
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"could not open file"

if not g.hasGrowthMatrix():
    raise IOError,"pickle file does not contain a growth matrix"

if not args.outfile is None:
    try:
        out = open(args.outfile,"w")
    except:
        out = sys.stdout
else:
    out = sys.stdout

if isinstance(g.growthmatrixgrid,int):
    for x in range(g.growthmatrixgrid):
        for y in range(g.growthmatrixgrid):
            print >> out,x,y,g.growthmatrix[x,y,0],g.growthmatrix[x,y,1]
        print >> out
else:
    raise NotImplementedError

if not out.closed:
    out.close()


