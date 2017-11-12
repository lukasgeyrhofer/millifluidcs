#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles/')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=False)

parser_pvd = parser.add_argument_group(description = "==== Parameters for interactions with PVD ====")
parser_pvd.add_argument("-Y","--PVD_Internal_Yield",type=float,nargs="*",default=[1,1])
parser_pvd.add_argument("-P","--PVD_Production",type=float,nargs="*",default=[1e-3,0])
parser_pvd.add_argument("-I","--PVD_Initial_Internal_Iron",type=float,nargs="*",default=[0,0])
parser_pvd.add_argument("-E","--PVD_Base_Iron_Influx",type=float,default=1)
parser_pvd.add_argument("-K","--PVD_Kpvd",type=float,default=1e-30)
parser_pvd.add_argument("-F","--PVD_Total_Iron",type=float,default=1e4)

parser_gm = parser.add_argument_group(description = "==== Growthmatrix calculation ====")
parser_gm.add_argument("-m","--maxM",default=40,type=int)
parser_gm.add_argument("-M","--stepM",default=1,type=int)

parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-o","--outfile",default=None)
parser_io.add_argument("-i","--infile",default=None)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")


args = parser.parse_args()

if args.outfile is None:
    raise IOError, "filename not specified"

if args.infile is None:
    g    = gc.GrowthDynamicsPyoverdin3(**vars(args))
    if args.verbose:print g.ParameterString()
    g.ComputeGrowthMatrix(size = args.maxM,step = args.stepM)

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





