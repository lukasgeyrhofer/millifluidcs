#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import inspect
import growthclasses as gc


def MakeGrowthDynList(module):
    gdl = dict()
    for x in inspect.getmembers(module):
        if x[0][:14] == 'GrowthDynamics' and inspect.isclass(x[1]):
            gdl[x[0]] = x[1]
    return gdl

def MakeDictFromParameters(params):
    p = dict()
    curkey = None
    curvalue = list()
    for entry in params:
        if is_number(entry):
            curvalue.append(float(entry))
        else:
            if not curkey is None:
                if len(curvalue) == 1:
                    p[curkey] = curvalue[0]
                elif len(curvalue) > 1:
                    p[curkey] = np.array(curvalue)
            curvalue = list()
            curkey = entry
    return p

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser_gc = parser.add_argument_group("==== GrowthDynamics ====")
parser_gc.add_argument("-d","--GrowthDynamics",default="")
parser_gc.add_argument("-L","--ParameterList",nargs="*",default=[])

parser_alg = parser.add_argument_group(description = "==== Parameters for algorithm ====")
parser_alg.add_argument("-t","--TimeIntegratorStep",default=1e-3,type=float)
parser_alg.add_argument("-O","--TimeIntegratorOutput",default=10,type=int)
parser_alg.add_argument("-M","--IntegrationMethod",choices = ['ownRK4','SciPy'],default = 'ownRK4')

parser_ic = parser.add_argument_group(description = "==== Initial conditions ====")
parser_ic.add_argument("-N","--initialconditions",default=[1],type=float,nargs="*")

args = parser.parse_args()

GrowthDynList = MakeGrowthDynList(gc)
param = MakeDictFromParameters(args.ParameterList)
param.update(**vars(args))

if 'GrowthDynamics' + args.GrowthDynamics in GrowthDynList.keys():
    g = GrowthDynList['GrowthDynamics' + args.GrowthDynamics](**param)
else:
    raise ValueError("'{}' not implemented in growthclasses.py".format('GrowthDynamics' + args.GrowthDynamics))


traj = g.Trajectory(args.initialconditions,TimeOutput=True)
for x in traj:
    print ' '.join(['{:14.6e}'.format(y) for y in x])

