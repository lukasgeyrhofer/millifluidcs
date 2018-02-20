#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import itertools


sys.path.append(sys.path[0] + '/..')
import growthclasses as gc


def correction_term(m1,m2,model,modelparameters):
    x = 0
    if m1+m2>0:x=float(m1)/(m1+m2)
    r = 1 + modelparameters['dy'][0] * (x-.5)
    if model == 'AB':
        if m1 * modelparameters['cmdline'][0] + m2 * modelparameters['cmdline'] < 1:
            r = 0
    elif model == 'PVD':
        r *= modelparameters['cmdline'][0]
    return r

def coord_to_inoc(c):
    return np.array([c[1]*c[0],(1-c[1])*c[0]])

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser_lattice = parser.add_argument_group(description = "==== Output lattice ====")
parser_lattice.add_argument("-N","--maxN",   default=50,  type=float)
parser_lattice.add_argument("-n","--stepN",  default=.1,  type=float)
parser_lattice.add_argument("-x","--stepsX", default=101, type=int)
parser_lattice.add_argument("-m","--maxM",   default=100, type=int)

parser_model = parser.add_argument_group(description = "==== Within droplet dynamics ====")
parser_model.add_argument("-M","--model",choices=['GY','PVD','AB'],default='GY')
parser_model.add_argument("-p","--modelparameters",nargs="*",type=float,default=[])

args = parser.parse_args()

g = gc.GrowthDynamics(**vars(args))

xlist = np.linspace( start = 0, stop = 1,         num  = args.stepsX)
nlist = np.arange(   start = 0, stop = args.maxN, step = args.stepN)
mlist = np.arange(   start = 0, stop = args.maxM,dtype=int)

y     = np.mean(g.yieldfactors)
dy    = (g.yieldfactors - y)/y
a     = np.mean(g.growthrates)
da    = (g.growthrates - a)/a

f1    = np.zeros((args.maxM,args.maxM))
f2    = np.zeros((args.maxM,args.maxM))
f3    = np.zeros((args.maxM,args.maxM))

modelparameters = {'cmdline':args.modelparameters,'dy':dy,'da':da}

for m1,m2 in itertools.product(mlist,repeat=2):
    n = float(m1+m2)
    if n > 0:
        x = m1/n
        f1[m1,m2] = x * np.power(n,-da[0]) * np.power(correction_term(m1,m2,args.model,modelparameters),1+da[0])
        f2[m1,m2] =     np.power(n, da[0]) * np.power(correction_term(m1,m2,args.model,modelparameters),1-da[0])
        f3[m1,m2] = x * np.power(n, da[0]) * np.power(correction_term(m1,m2,args.model,modelparameters),1-da[0])
    

lastx = -1
for coord in itertools.product(nlist,xlist):
    p = gc.PoissonSeedingVectors(mlist,coord_to_inoc(coord))

    avg_f1 = np.dot(p[1],np.dot(p[0],f1))
    avg_f2 = np.dot(p[1],np.dot(p[0],f2))
    avg_f3 = np.dot(p[1],np.dot(p[0],f2))
    
    if coord[1] < lastx:
        sys.stdout.write("\n")
    lastx = coord[1]

    sys.stdout.write("{:.3f} {:.6f} {:14.6e} {:14.6e} {:14.6e}\n".format(coord[0],coord[1],avg_f1,avg_f2,avg_f3))













