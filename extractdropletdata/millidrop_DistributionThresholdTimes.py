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

aparser = parser.add_argument_group(description = "==== algorithm parameters ====")
aparser.add_argument("-T","--threshold",type=float,default=.4)
aparser.add_argument("-C","--channel",type=str,default="Channel1_mean")

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-B","--bins",default=10,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)


args = parser.parse_args()

data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = args.splitBackForthTrajectories, datacolumns = ["time",args.channel], restrictionfile = args.restrictionfile)

for ExperimentLabel,Trajectories:
    thresholdtimes = list()
    for traj in Trajectories:
        t = traj[:,0]
        b = traj[:,1]
        
        thresholdtimes.append(np.min(t[b>args.threshold]))

    hist,bins = np.histogram(thresholdtimes,bins = args.bins,range = args.histogramrange,density = True)
    bins = bins[:-1] + np.diff(bins)/2.
    
    outfilename = args.outbasename + ExperimentLabel + ".threshold{:0.5e}".format(args.threshold)
    np.savetxt(outfilename,np.transpose(np.array([bins,hist])),fmt = "%.6e")


                    







