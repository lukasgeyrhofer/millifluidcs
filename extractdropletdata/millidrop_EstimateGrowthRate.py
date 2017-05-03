#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
ioparser = argparse.add_argument_group(description = "==== I/O parameters ====")
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-t","--templatefile",default=None)
parser.add_argument("-R"."--restrictionfile",default=None)

ffparser = parser.add_mutuall_exclusive_group(description = "==== Type of fit function ====")
ffparser.add_argument("-E","--exponential", default = False, action = "store_true")
ffparser.add_argument("-L","--logistic",    default = False, action = "store_true")

hparser argparse.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-B","--bins",default=10,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)

args = parser.parse_args()

data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = True)

if not args.restrictionfile is None:
    data.load_restrictions_from_file(args.restrictionfile)

growthrates = dict()

for experimentLabel, trajectories in data:
    if not growthrates.has_key(experimentLabel):
        growthrates[experimentLabel] = list()
    
    for trajectory in trajectories:
        t = trajectory[:,0] * 1e-3
        b = trajectory[:,1]
        
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
        r = (.99 * np.min(growthrates[experimentLabel]), 1.01 * np.max(growthrates[experimentLabel]))
    else:
        r = args.histogramrange
    h,b = np.histogram(growthrates[experimentLabel],bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    np.savetxt(experimentLabel + ".growthrates",np.transpose([b,h]))
