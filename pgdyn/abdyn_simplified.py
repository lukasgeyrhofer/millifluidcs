#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

def dyn(t,x,params):
    return np.array([
            params['lambda']*(1-x[1])*x[0],
            -params['eta']*x[1]*x[0]
           ])



parser = argparse.ArgumentParser()
parser.add_argument("-b","--antibiotics",default=2,type=float)
parser.add_argument("-n","--bacteria",default=10,type=float)
parser.add_argument("-l","--bactgrowth",default=1.5,type=float)
parser.add_argument("-e","--abdecay",default=1e-1,type=float)

parser.add_argument("-t","--integrationstep",type=float,default=1e-3)
parser.add_argument("-o","--outputstep",type=int,default=100)
parser.add_argument("-T","--maxtime",type=float,default=10)
args = parser.parse_args()


d = gc.TimeIntegrator(dynamics = dyn,initialconditions = np.array([args.bacteria,args.antibiotics]),step = args.integrationstep,params = {'lambda':args.bactgrowth,'eta':args.abdecay})
d.SetEndCondition("maxtime",args.maxtime)



print d.time,d
while not d.HasEnded():
    d.IntegrationStep(args.integrationstep * args.outputstep)
    print d.time,d


