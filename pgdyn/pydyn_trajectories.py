#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


def yieldfactor(pgconc):
    return params['yield'] + params['changeyield']/(1 + np.exp(-(pgconc - params['yieldhalfeffect'])/params['yieldspread'])


def dyn(t,xx,parameter):
    gr = params['growthrate']
    if xx[1] < 1e-200:
        gr = 0
    return np.array([   gr * xx[0],
                        -gr/yieldfactor(xx[1])*xx[0],
                        params['production']*xx[0] ])
                        
                       

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parserP = parser.add_argument_group(description = "other physiological parameters")
parserP.add_argument("-P","--production",type=float,default=1)
parserP.add_argument("-D","--degradation",type=float,default=1e-2)
parserP.add_argument("-d","--deltaY",type=float,default=1)
parserP.add_argument("-m","--yieldhalfeffect",type=float,default=1e3)
parserP.add_argument("-M","--yieldspread",type=float,default=2e2)
parserP.add_argument("-N","--initialN",type=float,default=10)
parserP.add_argument("-B","--initialB",type=float,default=.5)

parserA = parser.add_argument_group(description = "algorithm parameters")
parserA.add_argument("-s","--integrationstep",type=float,default=1e-3)
parserA.add_argument("-o","--outputstep",type=float,default=5e-2)

args = parser.parse_args()

global params
params = {  'growthrate'        : args.growthrates[0],
            'yield'             : args.yieldfactors[0],
            'yieldhalfeffect'   : args.yieldhalfeffect,
            'yieldspread'       : args.yieldspread,
            'production'        : args.production}

maxtime = args.mixingtime
if maxtime is None:
    maxtime = 20

x = np.zeros(3)
x[0] = args.initialN
x[1] = args.substrateconcentration

d = gc.TimeIntegrator(dynamics = dyn,initialconditions = x,step = args.integrationstep)
while d.time <= maxtime:
    print "{:.2f} {:11.4e} {:11.4e} {:11.4e}".format(d.time,*d.populations)
    d.IntegrationStep(args.outputstep)
print "{:.2f} {:11.4e} {:11.4e} {:11.4e}".format(d.time,*d.populations)





    
