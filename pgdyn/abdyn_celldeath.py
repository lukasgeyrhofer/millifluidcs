#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


def deathrate(abconc):
    bk = np.power(abconc,params['kappa'])
    return (params['gamma']+1)*bk/(bk + params['gamma'])


def dyn(t,xx,parameter):
    gr = params['growthrate']
    if xx[2] < 1e-200:
        gr = 0
    return np.array([   gr * ( 1 - deathrate(xx[4]))* xx[0],       # 0: alive cells
                        gr * deathrate(xx[4]) * xx[0],             # 1: dead cells
                        -gr * xx[0]/params['yield'],               # 2: substrate
                        params['production'] * xx[0],              # 3: production public good
                        -params['degradation'] * xx[3] * xx[4] ])  # 4: antibiotics
                       

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parserP = parser.add_argument_group(description = "other physiological parameters")
parserP.add_argument("-P","--production",type=float,default=1)
parserP.add_argument("-D","--degradation",type=float,default=1e-2)
parserP.add_argument("-K","--kappa",type=float,default=1)
parserP.add_argument("-G","--gamma",type=float,default=5)
parserP.add_argument("-N","--initialN",type=float,default=10)
parserP.add_argument("-B","--initialB",type=float,default=.5)

parserA = parser.add_argument_group(description = "algorithm parameters")
parserA.add_argument("-s","--integrationstep",type=float,default=1e-3)
parserA.add_argument("-o","--outputstep",type=float,default=5e-2)

args = parser.parse_args()

global params
params = {  'growthrate'    : args.growthrates[0],
            'yield'         : args.yieldfactors[0],
            'kappa'         : args.kappa,
            'gamma'         : args.gamma,
            'production'    : args.production,
            'degradation'   : args.degradation  }

maxtime = args.mixingtime
if maxtime is None:
    maxtime = 20

x = np.zeros(5)
x[0] = args.initialN
x[2] = args.substrateconcentration
x[4] = args.initialB

d = gc.TimeIntegrator(dynamics = dyn,initialconditions = x,step = args.integrationstep)
while d.time <= maxtime:
    print "{:.2f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}".format(d.time,*d.populations)
    d.IntegrationStep(args.outputstep)
    if d[0] < 1:d[0]=0
print "{:.2f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}".format(d.time,*d.populations)





    
