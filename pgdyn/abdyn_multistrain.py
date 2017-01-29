#!/usr/bin/env python


import argparse
import numpy as np
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

# 0     population producers
# 1     antibiotics inside producers
# 2     lactamase inside producers
# 3     population non-producers
# 4     antibiotics inside non-producers
# 5     antibiotics outside
# 6     lactamase outside
# 7     nutrients

def deathrate(rate,abconc,population):
    return rate * abconc * population
def degradation(rate,abconc,blconc):
    return rate * abconc * blconc
def production(rate,conc):
    return rate # only pure production?
def influx(rate,conc):
    return rate*conc
def outflux(rate,conc):
    return rate * conc


def dyn(t,x,params):
    return np.array([
        params['growth'] * x[0] - deathrate(params['deathrate'],x[1],x[0]),
        -params['growth'] * x[1] - degradation(params['degradation'],x[1],x[2]) + influx(params['influx'],x[5]) - outflux(params['outflux'],x[1]),
        -params['growth'] * x[2] + params['production'] - params['excretion'] * x[2],
        params['growth'] * x[3] - deathrate(params['deathrate'],x[4],x[3]),
        -params['growth'] * x[4] + influx(params['influx'],x[5]) - outflux(params['outflux'],x[4]),
        (- influx(params['influx'],x[5]) + outflux(params['outflux'],x[1])) * x[0] + (-influx(params['influx'],x[5]) + outflux(params['outflux'],x[4])) * x[3] - degradation(params['degradation'],x[5],x[6]),
        params['excretion'] * x[2] * x[0],
        -params['growth'] * x[0] / params['yield'] - params['growth'] * x[3] / params['yield']
    ])


parser = argparse.ArgumentParser()
parser.add_argument("-a","--growth",type = float,default = 1)
parser.add_argument("-y","--yield",type=float,default = 1)
parser.add_argument("-d","--deathrate",type=float,default = 1)
parser.add_argument("-r","--degradation",type=float,default=1)
parser.add_argument("-e","--excretion",type=float,default=1)
parser.add_argument("-i","--influx",type=float,default=1)
parser.add_argument("-o","--outflux",type=float,default=1)
parser.add_argument("-p","--production",type=float,default=1)

parser.add_argument("-S","--timestep",type=float,default=1e-3)
parser.add_argument("-M","--maxtime",type=float,default=20)
parser.add_argument("-N","--initialnonproducerbacteria",type=float,default=1e3)
parser.add_argument("-P","--initialproducerbacteria",type=float,default=1e3)
parser.add_argument("-B","--initialantibiotics",type=float,default=2)
parser.add_argument("-I","--initialnutrients",type=float,default=1e5)

args = parser.parse_args()

params = vars(args)
x = np.array([args.initialproducerbacteria,0,0,args.initialnonproducerbacteria,0,args.initialantibiotics,0,args.initialnutrients])
d = gc.TimeIntegrator(dynamics = dyn, initialconditions = x,step = args.timestep,params = params)
d.SetEndCondition("maxtime",args.maxtime)

while not d.HasEnded():
    print d
    
