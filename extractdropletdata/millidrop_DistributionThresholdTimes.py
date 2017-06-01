#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
ioparser.add_argument("-i","--infiles",nargs="*")
ioparser.add_argument("-t","--templatefile",default=None)
ioparser.add_argument("-r","--restrictionfile",default=None)
ioparser.add_argument("-o","--outbasename",default=None)
ioparser.add_argument("-B","--splitBackForthTrajectories",default=True,action="store_false")
ioparser.add_argument("-u","--timerescale",default=3.6e3,type=float)
ioparser.add_argument("-v","--verbose",default=False,action="store_true")

aparser = parser.add_argument_group(description = "==== algorithm parameters ====")
aparser.add_argument("-T","--threshold",type=float,default=.4)
aparser.add_argument("-C","--channel",type=str,default="Channel1_mean")

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-b","--bins",default=10,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)


args = parser.parse_args()

data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = args.splitBackForthTrajectories, datacolumns = ["time",args.channel], restrictionfile = args.restrictionfile)

for ExperimentLabel,Trajectories in data:
    thresholdtimes = list()
    for traj in Trajectories:
        t = traj[:,0] / args.timerescale
        b = traj[:,1]
        
        timesabovethreshold = t[b>=args.threshold]
        if len(timesabovethreshold) > 0:
            thresholdtimes.append(np.min(timesabovethreshold))
        
    if len(thresholdtimes) > 0:
        hist,bins = np.histogram(thresholdtimes,bins = args.bins,range = args.histogramrange,density = True)
        bins = bins[:-1] + np.diff(bins)/2.
        
        if args.verbose:
            print "{:15s}: storing histogram of {:d} values".format(ExperimentLabel,len(thresholdtimes))
        
        outfilename = args.outbasename + ExperimentLabel + ".threshold{:0.3e}".format(args.threshold)
        np.savetxt(outfilename,np.transpose(np.array([bins,hist])),fmt = "%.6e")


                    







