#!/usr/bin/env python

import numpy as np
import argparse
import math,sys

import millidrop_dataclass as mdc

from scipy.optimize import curve_fit

def logisticgrowth(t,a0,n0,k0,t0):
    return n0*np.exp(a0*(t-t0)) / (1.+n0/k0*(np.exp(a0*(t-t0))-1.))


parser = argparse.ArgumentParser()
parser.add_argument("-T","--templatefile")
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-l","--lowerthreshold", type=float, default = .027)
parser.add_argument("-M","--maxfev",type=int,default=1000)
args = parser.parse_args()


data = mdc.DropletData(templatefile = args.templatefile, infiles = args.infiles, splitBackForthTrajectories = True)


for label,traj in data:
    for t in traj:
        ic = np.array([.05,1e-3,4,10])
        
        fit0,fit1 = curve_fit(logisticgrowth,t[:,0]*1e3,t[:,1],p0=ic,maxfev = args.maxfev)
        print fit0
        
        exit(1)









