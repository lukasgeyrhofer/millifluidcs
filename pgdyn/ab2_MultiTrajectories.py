#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
parser_io.add_argument("-i","--infile")
parser_io.add_argument("-o","--outbasename",default=None)
parser_io.add_argument("-v","--verbose",default=False,action="store_true")
parser_io.add_argument("-t","--TimeIntegratorStep",type=float,default=1e-3)
parser_io.add_argument("-O","--OutputStep",type=float,default=1e-1)


parser_ab = parser.add_argument_group(description = "==== Dynamics of AB ====")
parser_ab.add_argument("-B","--ABconc",type=float,default=1.2)
parser_ab.add_argument("-g","--gamma",type=float,default=2)
parser_ab.add_argument("-k","--kappa",type=float,default=2)
parser_ab.add_argument("-p","--AB_Production_Efficiency",nargs="*",default=[1e-3,0])

parser_ic = parser.add_argument_group(description = "==== Initial conditions ====")
parser_ic.add_argument("-n","--minStrain1",type=float,default=10)
parser_ic.add_argument("-N","--maxStrain1",type=float,default=None)
parser_ic.add_argument("-z","--stepsStrain1",type=int,default=10)
parser_ic.add_argument("-l","--logStrain1",default=False,action="store_true")
parser_ic.add_argument("-m","--minStrain2",type=float,default=10)
parser_ic.add_argument("-M","--maxStrain2",type=float,default=None)
parser_ic.add_argument("-K","--stepsStrain2",type=int,default=10)
parser_ic.add_argument("-Z","--logStrain2",default=False,action="store_true")

args = parser.parse_args()

g = gc.GrowthDynamicsAntibiotics2(**vars(args))

if not args.maxStrain1 is None:
    if args.logStrain1: nlist1 = np.exp(np.linspace(start = np.log(args.minStrain1), stop = np.log(args.maxStrain1), num = args.stepsStrain1))
    else:               nlist1 = np.linspace(start = args.minStrain1, stop = args.maxStrain1, num = args.stepsStrain1)
else:                   nlist1 = np.array([args.minStrain1])
    
if not args.maxStrain2 is None:
    if args.logStrain2: nlist2 = np.exp(np.linspace(start = np.log(args.minStrain2), stop = np.log(args.maxStrain2), num = args.stepsStrain2))
    else:               nlist2 = np.linspace(start = args.minStrain2, stop = args.maxStrain2, num = args.stepsStrain2)
else:                   nlist2 = np.array([args.minStrain2])

for n1 in nlist1:
    for n2 in nlist2:
        if args.verbose:
            print n1,n2
        traj = g.Trajectory(args.OutputStep,initialconditions = np.array([n1,n2]))
        np.savetxt("{:s}_{:06.2f}_{:06.2f}".format(args.outbasename,n1,n2),traj)















