#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)


parser.add_argument("-S","--lowerSignalCutoff",default=None,type=float)
parser.add_argument("-R","--histogramRange",default = None,nargs=2,type=float)
parser.add_argument("-b","--histogramBins",default=None,type=int)
parser.add_argument("-v","--verbose",default = False,action="store_true")
args = parser.parse_args()

data = mdc.DropletData(**vars(args))

hRange = args.histogramRange
hBins  = args.histogramBins


for experimentlabel,trajectories in data:
    finalSignal = list()
    for trajectory in trajectories:
        if not args.lowerSignalCutoff is None:
            if trajectory[-1] > args.lowerSignalCutoff:
                finalSignal.append(trajectory[-1])
        else:
            finalSignal.append(trajectory[-1])
    
    outfilename = data.outbasename + "Yield_" + experimentlabel
    h,b = np.histogram(finalSignal,bins = hBins,range = hRange)
    b = b[:-1] + .5 * np.diff(b)
    np.savetxt(np.transpose([b,h]),outfilename)
    if args.verbose:
        





