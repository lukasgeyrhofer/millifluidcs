#!/usr/bin/env python3


import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc

parser = argparse.ArgumentParser()
parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
parser_io.add_argument("-o","--outfile",default=None,required = True)
parser_io.add_argument("-i","--infile",default=None)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")

parser = gc.AddGrowthDynamicsArguments(parser)
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=False)

parser_gm = parser.add_argument_group(description = "==== Parameters for growthmatrix ====")
parser_gm.add_argument("-m","--maxsize",type=int,default=100)
parser_gm.add_argument("-M","--step",type=int,default=1)

args = parser.parse_args()


if args.infile is None:
    g = gc.AssignGrowthDynamics(**vars(args))
    if args.verbose:print(g.ParameterString())
    g.ComputeGrowthMatrix(size = args.maxsize,step = args.step)

else:
    try:
        g = pickle.load(open(args.infile))
    except:
        raise IOError("could not open pickle file")
    if args.verbose:print g.ParameterString()
    if g.hasGrowthMatrix():
        g.ExtendGrowthMatrix(size = args.maxsize)
    else:
        raise IOError("pickle file does not contain growthmatrix")
    
try:
    fp = open(args.outfile,"w")
    pickle.dump(g,fp)
except:
    raise IOError("could not open file for pickle dump")
        
    


