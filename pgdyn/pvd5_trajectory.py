#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle
import inspect

sys.path.append(sys.path[0] + '/../mixingcycles/')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=False)

parser_pvd = parser.add_argument_group(description = "==== Parameters for interactions with PVD ====")
parser_pvd.add_argument("-Y","--PVD_Yield_Increase_Factor",type=float,default=2)
parser_pvd.add_argument("-P","--PVD_Production",type=float,nargs="*",default=[1e-3,0])

parser_ic = parser.add_argument_group(description = "==== Initital conditions ====")
parser_ic.add_argument("-N","--initialconditions",type=float,nargs="*",default=[1,1])
 
parser_io = parser.add_argument_group(description = "==== I/O ====")
#parser_io.add_argument("-o","--outfile",default=None,required = True)
#parser_io.add_argument("-i","--infile",default=None)
#parser_io.add_argument("-v","--verbose",default=False,action="store_true")
parser_io.add_argument("-t","--TimeIntegratorStep",type=float,default=1e-3)

args = parser.parse_args()

g = gc.GrowthDynamicsPyoverdin5(**vars(args))

traj = g.Trajectory(args.initialconditions,TimeOutput=True)

for x in traj:
    print ' '.join(['{:14.6f}'.format(y) for y in x])

