#!/usr/bin/env python

import numpy as np
import argparse
import sys


def get_growth(initialcells,parameters):
    # Estimate time when resources are used up. This is an important quantity when mixtures of strains coexist in droplet.
    # Usually the fast growth strain uses up a large share of the whole amount of resources
    def get_time_to_substrate_depletion(alpha=1.,maxsteps = 10000, precision = 1e-10):
        # internal function to have all parameters.
        t0 = 0
        if np.sum(initialcells) > 0:
            # initial condition for iteration is assuming only strain with highest expected yield is present
            p = (1.*initialcells/parameters['yield']).argmax() # get ID of strain
            t1 = np.log(parameters['substrate']*parameters['yield'][p]/initialcells[p]+1.)/parameters['growth'][p] # set first estimate of time for single strain
            i = 0
            while ((t1-t0)/t1)**2 > precision**2:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += alpha*(parameters['substrate']-np.sum(initialcells/parameters['yield']*(np.exp(parameters['growth']*t1)-1.)))/(np.sum(initialcells/parameters['yield']*parameters['growth']*np.exp(parameters['growth']*t1)))
                i+=1
                # should not iterate infinitely
                if i > maxsteps:
                    raise ValueError
            return min(t1,parameters['mixingtime'])
        else:
            return 0

    t = get_time_to_substrate_depletion()
    return initialcells * np.exp(parameters['growth'] * t) * parameters['dilution'],t



parser = argparse.ArgumentParser()
parser.add_argument("-m","--minM",type=int,default=1)
parser.add_argument("-M","--maxM",type=int,default=100)
parser.add_argument("-D","--stepM",type=int,default=1)
parser.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
parser.add_argument("-Y","--yieldfactor",type=float,nargs="*",default=[1.,2.])
parser.add_argument("-S","--substrateconcentration",type=float,default=1e4)
parser.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
parser.add_argument("-T","--mixingtime",type=float,default=12.)
args = parser.parse_args()



# generate parameter dictionary
p = {}
p['yield']      = np.array(args.yieldfactor)
p['growth']     = np.array(args.growthrates)
p['substrate']  = args.substrateconcentration
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilutionfactor

# get fixed point for single strain
a = p['growth'].argmax()
fp = np.zeros(len(p['growth']))
if p['mixingtime'] > 1/p['growth'][a]*np.log(1./p['dilution']):
    fp[a] = p['dilution']/(1.-p['dilution']) * p['yield'][a] * p['substrate']

for i in range(args.minM,args.maxM,args.stepM):
    c1 = np.array([i,1.])
    c0 = np.array([i,0.])
    n1,t1 = get_growth(c1,p)
    n0,t0 = get_growth(c0,p)
    print '{:5d} {:.12e} {:.12e} {:.12e} {:.12e} {:.12e}'.format(i,n1[0],n1[1],n0[0],n0[1],t1,t0)

