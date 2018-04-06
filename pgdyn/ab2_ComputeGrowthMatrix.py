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

parser_ab = parser.add_argument_group(description = "==== Dynamics of AB ====")
parser_ab.add_argument("-B","--ABconc",type=float,default=1.2)
parser_ab.add_argument("-g","--gamma",type=float,default=2)
parser_ab.add_argument("-k","--kappa",type=float,default=2)
parser_ab.add_argument("-p","--AB_Production_Efficiency",nargs="*",default=[1e-3,0])

parser_iterationmap = parser.add_argument_group(description = "==== Iterationmap ====")
parser_iterationmap.add_argument("-m","--maxsize",type=int,default=100)
parser_iterationmap.add_argument("-M","--step",type=int,default=1)

parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-o","--outfile",default=None)
parser_io.add_argument("-i","--infile",default=None)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")
parser_io.add_argument("-t","--TimeIntegratorStep",type=float,default=1e-3)
parser_io.add_argument("-O","--TimeIntegratorOutput",type=int,default=10)

args = parser.parse_args()

if args.outfile is None:
    raise IOError, "filename not specified"

if args.infile is None:
    g    = gc.GrowthDynamicsAntibiotics2(**vars(args))
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
        
    


