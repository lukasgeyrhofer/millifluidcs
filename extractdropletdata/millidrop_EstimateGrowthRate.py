#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math
from scipy.optimize import curve_fit

import millidrop_dataclass as mdc


def logisticgrowth(t,la0,ln0,lk0,t0):
    return np.exp(ln0 + np.exp(la0)*(t-t0)) / (1.+np.exp(ln0)*(np.exp(np.exp(la0)*(t-t0))-1.)/np.exp(lk0))


parser = argparse.ArgumentParser()
ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
ioparser.add_argument("-i","--infiles",nargs="*")
ioparser.add_argument("-t","--templatefile",default=None)
ioparser.add_argument("-r","--restrictionfile",default=None)


aparser = parser.add_argument_group(description = "==== Algorithm parameters ====")
aparser.add_argument("-m","--maxfev",default=5000,type=int)

ffparser = aparser.add_mutually_exclusive_group()
ffparser.add_argument("-E","--exponential", default = False, action = "store_true")
ffparser.add_argument("-L","--logistic",    default = False, action = "store_true")

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-B","--bins",default=10,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)

args = parser.parse_args()

data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = True)

if not args.restrictionfile is None:
    data.load_restrictions_from_file(args.restrictionfile)

if not args.exponential and not args.logistic:
    mode = "exponential"
else:
    mode = "exponential" if args.exponential else "logistic"

growthrates = dict()

for experimentLabel, trajectories in data:
    if not growthrates.has_key(experimentLabel):
        growthrates[experimentLabel] = list()
    
    for trajectory in trajectories:
        t = trajectory[:,0] / 3600.
        b = trajectory[:,1]
        
        if mode == "exponential":
            if len(t) >= 2:
                sx  = np.sum(t)
                sxx = np.sum(t*t)
                sy  = np.sum(np.log(b))
                sxy = np.sum(t * np.log(b))
                n   = len(t)
                gr  = (n * sxy - sx * sy)/(n*sxx - sx*sx)

                growthrates[experimentLabel].append(gr)
                
        elif mode == "logistic":
            # length of trajectory has to contain more points than number of fitting parameters
            if len(t) > 4:
                ic = np.array([np.log(abs(np.log(b[0]/b[1])/(t[0] - t[1]))),np.log(b[0]),np.log(b[-1]),t[0]])
                fitMEAN,fitCOV = curve_fit(logisticgrowth,t,b,p0 = ic,maxfev = args.maxfev)
                growthrates[experimentLabel].append(np.exp(fitMEAN[0]))

    
    print "{:15s} {:.4f} (± {:.4f}) 1/h".format(experimentLabel,np.mean(growthrates[experimentLabel]),np.sqrt(np.std(growthrates[experimentLabel])))

    if args.histogramrange is None:
        r = (.99 * np.min(growthrates[experimentLabel]), 1.01 * np.max(growthrates[experimentLabel]))
    else:
        r = args.histogramrange
    h,b = np.histogram(growthrates[experimentLabel],bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    np.savetxt(experimentLabel + ".growthrates",np.transpose([b,h]))
