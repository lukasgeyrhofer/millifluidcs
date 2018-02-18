#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


def beta(abconc):
    bk = np.power(abconc,args.kappa)
    return 1 - (1+args.logkill)*bk/(bk + args.logkill)

def monod(substr):
    return substr/(args.monodKC + substr)


def dynAB(t,x,params):
    a = args.growthrates  * beta(x[-1]) # * monod(x[-3])
    if x[-3] == 0:
        a = np.zeros(len(args.growthrates))
    if np.any(x[:-3] < 1):
        (x[:-3])[x[:-3] < 1] = 0
    return np.concatenate([ a*x[:-3],                                          # growth of strains
                            np.array( [ -np.sum(a/args.yieldfactors*x[:-3]),   # decay of nutrients
                                        np.sum(args.PGproduction*x[:-3]),      # production of public good
                                        -args.PGreductionAB*x[-1]*x[-2] ])])   # reduction of antibiotics by public good


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime = 24)

parserAB = parser.add_argument_group(description = "Parameters for interactions with antibiotics")
parserAB.add_argument("-k","--kappa",type=float,default=1)
parserAB.add_argument("-l","--logkill",type=float,default=2)
parserAB.add_argument("-P","--PGproduction",nargs="*",default=[1,0])
parserAB.add_argument("-R","--PGreductionAB",type=float,default=1)
parserAB.add_argument("-K","--monodKC",type=float,default=1e2)

parserALG = parser.add_argument_group(description = "Parameters for algorithm")
parserALG.add_argument("-t","--integrationstep",default=1e-3,type=float)
parserALG.add_argument("-o","--outputstep",default=100,type=int)

parserIC = parser.add_argument_group(description = "Initial conditions")
parserIC.add_argument("-B","--ABconc",type=float,default=.5)
parserIC.add_argument("-N","--initialN",default=1,type=float)

global args
args = parser.parse_args()

args.growthrates = np.array(args.growthrates,dtype=np.float64)
args.yieldfactors  = np.array(args.yieldfactors,dtype=np.float64)
args.PGproduction = np.array(args.PGproduction,dtype=np.float64)
args.PGreductionAB = np.array(args.PGreductionAB,dtype=np.float64)
numstrains = len(args.growthrates)

ic = np.concatenate([np.ones(numstrains) * args.initialN, np.array([args.substrateconcentration,0,args.ABconc])])

d = gc.TimeIntegrator(  dynamics = dynAB ,
                        initialconditions = ic,
                        step = args.integrationstep)
#d.SetEndCondition("reachzero",numstrains)
d.SetEndCondition("maxtime",args.mixingtime)

print "%6.3f"%d.time,d
while not d.HasEnded():
    d.IntegrationStep(args.integrationstep * args.outputstep)
    print "%6.3f"%d.time,d






