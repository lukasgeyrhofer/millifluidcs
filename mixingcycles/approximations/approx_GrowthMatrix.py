#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/..')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-o","--outfile",                        default = "out.pickle")
parser_io.add_argument("-i","--infile",                         default = None)
parser_io.add_argument("-v","--verbose", action = "store_true", default = False)

parser_iterationmap = parser.add_argument_group(description = "Parameters for iterationmap")
parser_iterationmap.add_argument("-m","--maxsize",type=int,default=100)
parser_iterationmap.add_argument("-k","--step",   type=int,default=1)

parser_model = parser.add_argument_group(description = "==== Within droplet dynamics ====")
parser_model.add_argument("-M","--model",           choices = ['GY','PVD','AB'], default = 'GY')
parser_model.add_argument("-p","--modelparameters", nargs = "*", type = float,   default = [])

args = parser.parse_args()


if args.outfile is None:
    raise IOError, "filename not specified"

if args.infile is None:
    g    = gc.GrowthDynamicsApprox(**vars(args))
    if args.verbose:print g.ParameterString()
    g.ComputeGrowthMatrix(size = args.maxsize,step = args.step)

else:
    try:
        g = pickle.load(open(args.infile))
    except:
        raise IOError,"could not open pickle file"
    if args.verbose:print g.ParameterString()
    if g.hasGrowthMatrix():
        g.ExtendGrowthMatrix(size = args.maxsize)
    else:
        raise IOError,"pickle file does not contain growthmatrix"
    
try:
    fp = open(args.outfile,"w")
    pickle.dump(g,fp)
except IOError:
    print >> sys.stderr,"could not open file for pickle dump"
    exit(1)

