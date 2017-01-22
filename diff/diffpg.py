#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

def secondder(p):
    pp = np.concatenate([p[1:],np.array([2*p[-1] - p[-2]])])
    pm = np.concatenate([np.array([2*p[0] - p[1]]),p[:-1]])
    return (pp - 2.*p + pm)/(dx*dx)
    

def dyn(t,p,params):
    return growth + secondder(p)


parser = argparse.ArgumentParser()
parser.add_argument("-a","--productionrate",type=float,default=1)
parser.add_argument("-D","--diffusionconstant",type=float,default=1)
parser.add_argument("-d","--dx",type=float,default=.1)
parser.add_argument("-z","--space0",type=int,default=250)
parser.add_argument("-s","--space",type=int,default=500)
parser.add_argument("-t","--integrationstep",type=float,default=1e-3)
parser.add_argument("-o","--outputstep",type=int,default=100)
parser.add_argument("-T","--maxtime",type=float,default=10)
args = parser.parse_args()



x = (np.arange(args.space)-args.space0)*args.dx

space  = args.space
space0 = args.space0
global dx
dx     = args.dx

p = np.zeros(args.space)

global growth
growth = np.zeros(args.space)
growth[space0] = args.productionrate


d = gc.TimeIntegrator(  dynamics = dyn,
                        initialconditions = np.zeros(space),
                        step = args.integrationstep)
d.SetEndCondition("maxtime",args.maxtime)

print "%6.3f"%d.time,d
while not d.HasEnded():
    d.IntegrationStep(args.integrationstep * args.outputstep)
    for a1,a2 in zip(x,d.populations):
        print "%6.3f"%d.time,a1,a2
    print




