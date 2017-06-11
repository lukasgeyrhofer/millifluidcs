#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


def dynPVD(t,x,params):
    p = params['PVDincreaseS'] if x[-1] <= params['substrate'] * params['PVDmaxFactorS'] else 0
    if x[-3] > 0:
        a = params['growthrates']
    else:
        a = np.zeros(len(x)-3)
    return np.concatenate([ a*x[:-3],   np.array([  np.sum(-a*x[:-3]/params['yieldfactors']) + p * x[-2],
                                                    np.sum(params['PVDproduction']*x[:-3]),
                                                    p * x[-2]   ])])


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime = 24)

parserPVD = parser.add_argument_group(description = "Parameters for dynamics of PVD")
parserPVD.add_argument("-k","--PVDincreaseS",type=float,default=.2)
parserPVD.add_argument("-K","--PVDmaxFactorS",type=float,default=1.2)
parserPVD.add_argument("-p","--PVDproduction",nargs="*",default=[1,0])

parserALG = parser.add_argument_group(description = "Parameters for algorithm")
parserALG.add_argument("-t","--integrationstep",default=1e-3,type=float)
parserALG.add_argument("-o","--outputstep",default=100,type=int)

parserIC = parser.add_argument_group(description = "Initial conditions")
parserIC.add_argument("-N","--initialN",default=1,type=float)
parserIC.add_argument("-V","--PVDconc",type=float,default=0.0)

args = parser.parse_args()
kwargs = vars(args)

PVDparams = {  'PVDproduction' :  np.array(kwargs.get("PVDproduction",np.zeros(len(args.growthrates))),dtype=np.float64),
               'PVDincreaseS'  :  kwargs.get("PVDincreaseS",1),
               'PVDmaxFactorS' :  kwargs.get("PVDmaxFactorS",1),   # initial condition PG concentration
               'PVDconc' :        kwargs.get("PVDconc",0)}  # initial concentration antibiotics measured in zMIC

PVDparams['growthrates'] = np.array(args.growthrates)
PVDparams['yieldfactors'] = np.array(args.yieldfactors)
PVDparams['substrate'] = args.substrateconcentration


ic = np.concatenate([np.repeat(args.initialN,len(args.growthrates)),np.array([args.substrateconcentration,PVDparams['PVDconc'],0])])

d = gc.TimeIntegrator(dynamics = dynPVD,initialconditions = ic,step = args.integrationstep,params = PVDparams)
d.SetEndCondition("maxtime",args.mixingtime)

while not d.HasEnded():
    d.IntegrationStep(args.integrationstep * args.outputstep)
    print "%6.3f"%d.time,d



