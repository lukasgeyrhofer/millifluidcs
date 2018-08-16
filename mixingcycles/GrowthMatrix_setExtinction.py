#!/usr/bin/env python3

import argparse
import pickle

import growthclasses as gc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile")
parser.add_argument("-T","--extinctionThreshold",type=float,default=1)
args = parser.parse_args()

g = pickle.load(open(args.infile,'rb'),encoding='bytes')
g.setGrowthMatrixValues(args.extinctionThreshold,0,'below')
pickle.dump(g,open(args.outfile,'wb'))
