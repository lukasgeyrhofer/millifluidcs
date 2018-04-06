#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime = 24)

parser_ab = parser.add_argument_group(description = "==== Dynamics of AB ====")
parser_ab.add_argument("-B","--ABconc",type=float,default=1.2)
parser_ab.add_argument("-g","--gamma",type=float,default=2)
parser_ab.add_argument("-k","--kappa",type=float,default=2)
parser_ab.add_argument("-p","--AB_Production_Efficiency",nargs="*",default=[1e-3,0])

parser_alg = parser.add_argument_group(description = "==== Parameters for algorithm ====")
parser_alg.add_argument("-t","--TimeIntegratorStep",default=1e-3,type=float)
parser_alg.add_argument("-O","--TimeIntegratorOutput",default=10,type=int)

parser_ic = parser.add_argument_group(description = "==== Initial conditions ====")
parser_ic.add_argument("-N","--initialconditions",default=[1],type=float,nargs="*")

args = parser.parse_args()


g = gc.GrowthDynamicsAntibiotics2(**vars(args))

traj = g.Trajectory(args.initialconditions,TimeOutput=True)

for x in traj:
    print ' '.join(['{:14.6e}'.format(y) for y in x])




