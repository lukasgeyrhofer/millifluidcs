#!/usr/bin/env python

import argparse
import numpy as np
import sys,math


sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

def abgrowth(b,gamma,kappa):
    bk = np.power(b,kappa)
    return 1.-(gamma+1.)*bk/(bk+gamma)


def dyn_env(t,x,params):
    return np.array([
        abgrowth(x[1],params['gamma'],params['kappa'])*x[0],
        -params['tau']*x[1]
        ])

def dyn_int(t,x,params):
    return np.array([
        abgrowth(x[1],params['gamma'],params['kappa'])*x[0],
        -params['tau']*x[1]*x[0]
        ])

def dyn_enz(t,x,params):
    return np.array([
        abgrowth(x[1],params['gamma'],params['kappa'])*x[0],
        -params['tau']*x[1]*x[2],
        params['tau']*x[0]
        ])

def dyn_bnd(t,x,params):
    return np.array([
        abgrowth(x[1],params['gamma'],params['kappa'])*x[0],
        params['tau']*(abgrowth(x[1],params['gamma'],params['kappa'])-1)*x[0]
        ])


parser = argparse.ArgumentParser()
parser.add_argument("-N","--initialbacteria",type=float,default=1e3)
parser.add_argument("-B","--initialantibiotics",type=float,default=1)

parser.add_argument("-g","--gamma",type=float,default=2)
parser.add_argument("-k","--kappa",type=float,default=2)
parser.add_argument("-t","--tau",type=float,default=1e-3)

parser.add_argument("-S","--integrationstep",type=float,default=1e-3)
parser.add_argument("-O","--outputstep",type=int,default=100)
parser.add_argument("-T","--maxtime",type=float,default=10)
parser.add_argument("-F","--outputbasefilename",default="out")
args = parser.parse_args()

d1 = gc.TimeIntegrator(dynamics = dyn_env,initialconditions = np.array([args.initialbacteria,args.initialantibiotics])  ,step = args.integrationstep,params = {'kappa':args.kappa,'gamma':args.gamma,'tau':args.tau})
d2 = gc.TimeIntegrator(dynamics = dyn_int,initialconditions = np.array([args.initialbacteria,args.initialantibiotics])  ,step = args.integrationstep,params = {'kappa':args.kappa,'gamma':args.gamma,'tau':args.tau})
d3 = gc.TimeIntegrator(dynamics = dyn_enz,initialconditions = np.array([args.initialbacteria,args.initialantibiotics,0]),step = args.integrationstep,params = {'kappa':args.kappa,'gamma':args.gamma,'tau':args.tau})
d4 = gc.TimeIntegrator(dynamics = dyn_bnd,initialconditions = np.array([args.initialbacteria,args.initialantibiotics])  ,step = args.integrationstep,params = {'kappa':args.kappa,'gamma':args.gamma,'tau':args.tau})

d1.SetEndCondition("maxtime",args.maxtime)
d2.SetEndCondition("maxtime",args.maxtime)
d3.SetEndCondition("maxtime",args.maxtime)
d4.SetEndCondition("maxtime",args.maxtime)

f1 = args.outputbasefilename + '_env'
fp = open(f1,'w')
print >> fp,d1.time,d1
while not d1.HasEnded():
    d1.IntegrationStep(args.integrationstep * args.outputstep)
    print >>fp, d1.time,d1
fp.close()

f2 = args.outputbasefilename + '_int'
fp = open(f2,'w')
print >> fp,d2.time,d2
while not d2.HasEnded():
    d2.IntegrationStep(args.integrationstep * args.outputstep)
    print >>fp, d2.time,d2
fp.close()

f3 = args.outputbasefilename + '_enz'
fp = open(f3,'w')
print >> fp,d3.time,d3
while not d3.HasEnded():
    d3.IntegrationStep(args.integrationstep * args.outputstep)
    print >>fp, d3.time,d3
fp.close()

f4 = args.outputbasefilename + '_bnd'
fp = open(f4,'w')
print >> fp,d4.time,d4
while not d4.HasEnded():
    d4.IntegrationStep(args.integrationstep * args.outputstep)
    print >>fp, d4.time,d4
fp.close()


