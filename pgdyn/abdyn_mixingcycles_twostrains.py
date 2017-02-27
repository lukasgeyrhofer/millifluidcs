#!/usr/bin/env python


import argparse
import numpy as np
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


def death(abconc,params):
    return (params['gamma'] + 1)*np.power(abconc,params['kappa'])/(np.power(abconc,params['kappa']) + params['gamma'])
def growthrate(substrate,params):
    if substrate < 1e-100:  return 0
    else:                   return params['growthrate']
def growth(abconc,substrate,params):
    return growthrate(substrate,params) * (1 - death(abconc,params))
    
def dyn(t,x,params):
    if x[0] < 1:    x[0] = 0
    if x[1] < 1:    x[1] = 0 
    return np.array([
        growth(x[3],x[2],params) * x[0],
        params['phi'] * growth(x[3],x[2],params) * x[1],
        -growthrate(x[2],params)*(x[0] + params['phi']*x[1])/params['yieldfactor'],
        -params['sigma'] * x[3] * x[1],
        growthrate[x[2],params] * death(x[3],params) * (x[0] + params['phi']*x[1])
    ])

parser = argparse.ArgumentParser()

parser_params = parser.add_argument_group(description = "=== Parameters of microbial interactions ===")
parser_params.add_argument("-a","--growthrate",type=float,default=1)
parser_params.add_argument("-y","--yieldfactor",type=float,default=1)
parser_params.add_argument("-g","--gamma",type=float,default=2)
parser_params.add_argument("-k","--kappa",type=float,default=2)
parser_params.add_argument("-s","--sigma",type=float,default=1e-3)
parser_params.add_argument("-p","--phi",type=float,default=1-1e-2)
    
parser_ic = parser.add_argument(description = "=== Initial conditions ===")
parser_ic.add_argument("-N","--populationsize",type=float,default=1e2)
parser_ic.add_argument("-B","--antibioticconcentration",type=float,default=2)
parser_ic.add_argument("-S","--substrate",type=float,default=1e5)
parser_ic.add_argument("-P","--fractionreducers",type=float,default=.5)

parser_alg = parser.add_argument_group(description = "=== Algorithm parameters ===")
parser_alg.add_argument("-T","--mixingtime",type=float,default=24.)
parser_alg.add_argument("-I","--integrationstep",type=float,default=1e-3)
parser_alg.add_argument("-O","--outputstep",type=int,default=100)

args = parser.parse_args()

d = gc.TimeIntegrator(dynamics = dyn,initialconditions = np.array([args.populationsize * (1-args.fractionreducers),args.populationsize*args.fractionreducers,args.substrate,args.antibioticconcentration,0]), params = {'growthrate':args.growthrate,'yieldfactor':args.yieldfactor,'gamma':args.gamma,'kappa':args.kappa,'sigma':args.sigma,'phi':args.phi},step = args.integrationstep)
d.SetEndCondition('maxtime',args.maxtime)

print d.time,d
while not d.HasEnded():
    d.IntegrationStep(args.integrationstep * args.outputstep)
    print d.time,d



