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

parserALG = parser.add_argument_group(description = "==== Parameters for algorithm ====")
parserALG.add_argument("-t","--integrationstep",default=1e-3,type=float)
parserALG.add_argument("-o","--outputstep",default=100,type=int)

parserIC = parser.add_argument_group(description = "==== Initial conditions ====")
parserIC.add_argument("-N","--initialN",default=[1],type=float,nargs="*")

args = parser.parse_args()


g = gc.GrowthDynamicsAntibiotics2(**vars(args))

ns  = g.numstrains
lic = len(args.initialN)
if lic < ns:    ic = np.concatenate([args.initialN,np.zeros(ns-lic)])
else:           ic = np.array(args.initialN[:ns])

t = g.Trajectory(args.integrationstep * args.outputstep,ic)

for x in t:
    print '{:.2f} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(*x)










