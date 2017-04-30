#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-T","--templatefile",default=None)
parser.add_argument("-m","--lowerthreshold",default=.02,type=float)
parser.add_argument("-M","--upperthreshold",default=2,type=float)
parser.add_argument("-B","--bins",default=10,type=int)
parser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)
args = parser.parse_args()

data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = True)

growthrates = dict()

for experimentLabel, trajectories in data:
    if not growthrates.has_key(experimentLabel):
        growthrates[experimentLabel] = list()
    
    for trajectory in trajectories:
        t = trajectory[:,0] * 1e-3
        b = trajectory[:,1]
        
        t = t[args.lowerthreshold < b]
        b = b[args.lowerthreshold < b]
        
        t = t[args.upperthreshold > b]
        b = b[args.upperthreshold > b]
        
        if len(t) >= 2:
            sx  = np.sum(t)
            sxx = np.sum(t*t)
            sy  = np.sum(np.log(b))
            sxy = np.sum(t * np.log(b))
            n   = len(t)
            gr  = (n * sxy - sx * sy)/(n*sxx - sx*sx)

            growthrates[experimentLabel].append(gr)
    
    print "{:15s} {:.4f} (± {:.4f}) 1/h".format(experimentLabel,np.mean(growthrates[experimentLabel]),np.sqrt(np.std(growthrates[experimentLabel])))

    if args.histogramrange is None:
        r = (.2,.7)
    else:
        r = args.histogramrange
    h,b = np.histogram(growthrates[experimentLabel],bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    np.savetxt(experimentLabel + ".growthrates",np.transpose([b,h]))
