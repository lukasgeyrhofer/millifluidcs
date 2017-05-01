#!/usr/bin/env python

import numpy as np
import argparse
import math,sys

import millidrop_dataclass as mdc

from scipy.optimize import curve_fit

def logisticgrowth(t,la0,ln0,lk0,t0):
    return np.exp(ln0 + np.exp(la0)*(t-t0)) / (1.+np.exp(ln0)*(np.exp(np.exp(la0)*(t-t0))-1.)/np.exp(lk0))


parser = argparse.ArgumentParser()
parser.add_argument("-T","--templatefile")
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-l","--lowerthreshold", type=float, default = .027)
parser.add_argument("-t","--maxtime", type = float, default= 140)
parser.add_argument("-M","--maxfev",type=int,default=5000)
parser.add_argument("-B","--bins",type=int,default=10)
parser.add_argument("-v","--verbose",default=False,action="store_true")
args = parser.parse_args()


data = mdc.DropletData(templatefile = args.templatefile, infiles = args.infiles, splitBackForthTrajectories = True)


growthrates = dict()
yields      = dict()

for label,trajectories in data:
    growthrates[label] = list()
    yields[label]      = list()
    i                  = 0
    for traj in trajectories:
        i += 1
        t  = traj[:,0] * 1e-3
        b  = traj[:,1]
        
        t  = t[b > args.lowerthreshold]
        b  = b[b > args.lowerthreshold]
        
        b = b[t < args.maxtime]
        t = t[t < args.maxtime]
        
        try:
            ic = np.array([.05,np.log(b[0]),np.log(b[-1]),t[0]])
            fit0,fit1 = curve_fit(logisticgrowth,t,b,p0=ic,maxfev = args.maxfev)
            #print fit0,ic,fit0 == ic
            
            growthrates[label].append(np.exp(fit0[0]))
            yields[label].append(np.exp(fit0[2]))
        
        except:
            if args.verbose:
                print "error with '{}', trajectory {}".format(label,i)
            continue
        
        
    
    r = [-.05,.25]
    h,b = np.histogram(growthrates[label],bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    np.savetxt(label + ".growthratesE",np.transpose([b,h]))

    r = [0,6]
    h,b = np.histogram(yields[label],bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    np.savetxt(label + ".yieldsE",np.transpose([b,h]))







