#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import millidrop_dataclass as mdc
parser = argparse.ArgumentParser()
ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
ioparser.add_argument("-i","--infiles",nargs="*")
ioparser.add_argument("-t","--templatefile",default=None)
ioparser.add_argument("-r","--restrictionfile",default=None)
ioparser.add_argument("-o","--outbasename",default='')
ioparser.add_argument("-V","--write_values_to_outfile",default=False,action = "store_true")
ioparser.add_argument("-v","--verbose",default=False,action="store_true")

aparser = parser.add_argument_group(description = "==== Algorithm parameters ====")

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-B","--bins",default=120,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=[-.1,1.1])

args = parser.parse_args()
data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = True)





for experimentLabel, trajectories in data:
    instantgrowthrate = list()
    for trajectory in trajectories:
        for i in range(1,len(trajectory)):
            instantgrowthrate.append(np.log(trajectory[i,1]/trajectory[i-1,1])/(trajectory[i,0]-trajectory[i-1,0])*1e3)
    instantgrowthrate = np.array(instantgrowthrate)
    if args.histogramrange is None:
        r = (.99 * np.min(instantgrowthrate), 1.01 * np.max(instantgrowthrate))
    else:
        r = args.histogramrange
    h,b = np.histogram(instantgrowthrate,bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    outfilename = args.outbasename + experimentLabel + ".instantgrowthrates"
    np.savetxt(outfilename,np.transpose([b,h]),fmt = '%.6e')
    

