#!/usr/bin/env python

import argparse
import pickle

import growthclasses as gc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile")
parser.add_argument("-T","--extinctionThreshold",type=float,default=1)
args = parser.parse_args()

g = pickle.load(open(args.infile))
g.setExtinctionThreshold(args.extinctionThreshold)
pickle.dump(g,open(args.outfile,"w"))
