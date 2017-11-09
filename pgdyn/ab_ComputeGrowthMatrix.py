#!/usr/bin/env python


import numpy as np
import argparse
import sys,math
import pickle
import inspect

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=False)

parser_ab = parser.add_argument_group(description = "Parameters for dynamics of AB")
parser_ab.add_argument("-B","--ABconc",type=float,default=1.2)
parser_ab.add_argument("-g","--gamma",type=float,default=2)
parser_ab.add_argument("-k","--kappa",type=float,default=2)
parser_ab.add_argument("-p","--PGproduction",nargs="*",default=[1,0])
parser_ab.add_argument("-r","--PGreductionAB",type=float,default=1e-3)

parser_iterationmap = parser.add_argument_group(description = "Parameters for iterationmap")
parser_iterationmap.add_argument("-m","--maxsize",type=int,default=100)
parser_iterationmap.add_argument("-M","--step",type=int,default=1)

parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-i","--infile",default=None)
parser.add_argument("-v","--verbose",default=False,action="store_true")

args = parser.parse_args()

if args.outfile is None:
    raise IOError, "filename not specified"

if args.infile is None:
    g    = gc.GrowthDynamicsAntibiotics(**vars(args))
    if args.verbose:print g.ParameterString()
    g.ComputeGrowthMatrix(size = args.maxsize,step = args.step)

else:
    try:
        g = pickle.load(open(args.infile))
    except:
        raise IOError,"could not open pickle file"
    if args.verbose:print g.ParameterString()
    if g.hasGrowthMatrix():
        g.ExtendGrowthMatrix(size = args.maxsize,step = args.step)
    else:
        raise IOError,"pickle file does not contain growthmatrix"
    
try:
    fp = open(args.outfile,"w")
    pickle.dump(g,fp)
    fp.close()
except IOError:
    raise IOError("could not open file for pickle dump")
        
    


