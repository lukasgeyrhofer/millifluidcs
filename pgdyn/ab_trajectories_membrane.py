#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


def delta(abconc):
    bk = np.power(abconc,args.kappa)
    return args.growthrates * (1+args.logkill)*bk/(bk + args.logkill)


def dynAB(t,x,params):
    if np.any(x[:numstrains] < 1):
        (x[:numstrains])[x[:numstrains] < 1] = 0
    growthrates = args.growthrates
    if x[-1] <= 0:
        growthrates = np.zeros(numstrains)
        
        
        
    return np.concatenate([                                                                                                                                                                # growth of strains
                (growthrates - delta(x[numstrains:2*numstrains])) * x[:numstrains],                                                                                                        # ab, inside membrane
                args.ABdiffusion[0] * x[-3] - (args.ABdiffusion[1] + growthrates + args.ABreduction * x[2*numstrains:3*numstrains]) * x[numstrains:2*numstrains],                          # pg, inside membrane
                args.PGdiffusion[0] * x[-2] - (args.PGdiffusion[1] + growthrates) * x[2*numstrains:3*numstrains] + args.PGproduction,
                np.array([                                                                                                                                                                 # ab, outside membrane
                    -args.ABdiffusion[0] * x[-3] * np.sum(x[:numstrains]) + args.ABdiffusion[1] * np.dot(x[:numstrains], x[  numstrains:2*numstrains]) - args.ABreduction * x[-3] * x[-2], # pg, outside membrane
                    -args.PGdiffusion[0] * x[-2] * np.sum(x[:numstrains]) + args.PGdiffusion[1] * np.dot(x[:numstrains], x[2*numstrains:3*numstrains]),                                    # substrate
                    -np.sum(growthrates * x[:numstrains]/args.yieldfactors)
                    ])
                ])


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime = 24)

parserAB = parser.add_argument_group(description = "Parameters for interactions with antibiotics")
parserAB.add_argument("-k","--kappa",type=float,default=1)
parserAB.add_argument("-l","--logkill",type=float,default=2)
parserAB.add_argument("-B","--ABconc0",type=float,default=.5)
parserAB.add_argument("-R","--ABreduction",type=float,default=1)
parserAB.add_argument("-D","--ABdiffusion",type=float,nargs=2,default=[1,1])
parserAB.add_argument("-P","--PGproduction",nargs="*",default=[1,0])
parserAB.add_argument("-d","--PGdiffusion",type=float,nargs=2,default=[1,1])
parserAB.add_argument("-K","--monodKC",type=float,default=1e2)

parser.add_argument("-t","--integrationstep",default=1e-3,type=float)
parser.add_argument("-o","--outputstep",default=100,type=int)
parser.add_argument("-N","--initialN",default=1,type=float)
#parser.add_argument("-B","--initialB",default=.5,type=float)

global args
args = parser.parse_args()

args.growthrates   = np.array(args.growthrates,  dtype=np.float64)
args.yieldfactors  = np.array(args.yieldfactors, dtype=np.float64)
args.PGproduction  = np.array(args.PGproduction, dtype=np.float64)

global numstrains
numstrains = len(args.growthrates)

ic = np.concatenate([np.ones(numstrains) * args.initialN, np.zeros(2*numstrains), np.array([args.ABconc0,0,args.substrateconcentration])])

d = gc.TimeIntegrator(  dynamics = dynAB ,
                        initialconditions = ic,
                        step = args.integrationstep)
#d.SetEndCondition("reachzero",numstrains)
d.SetEndCondition("maxtime",args.mixingtime)

print "%6.3f"%d.time,d
while not d.HasEnded():
    d.IntegrationStep(args.integrationstep * args.outputstep)
    print "%6.3f"%d.time,d






