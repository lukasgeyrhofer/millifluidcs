#!/usr/bin/env python3

import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/..')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-O","--OverwriteInPlace",default=False,action="store_true")
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile,'rb'),encoding = 'bytes')
except:
    raise IOError("could not open pickle file '{}'".format(args.infile))

if args.outfile is None:
    if args.OverwriteInPlace:
        # overwrite file in place
        fp = open(args.infile,'wb')
    else:
        raise IOError("no output file specified")
else:
    fp = open(args.outfile,'wb')

pickle.dump(g,fp)
fp.close()

