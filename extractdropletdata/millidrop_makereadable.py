#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-t","--templatefile",default=None)
parser.add_argument("-C","--columns",nargs="*",default=["Channel1_mean"])
parser.add_argument("-B","--splitBackForthTrajectories",default=False,action="store_true")
args = parser.parse_args()

if not "time" in args.columns:
    columns = ["time"] + args.columns
else:
    columns = args.columns

data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, datacolumns = columns, splitBackForthTrajectories = True)

for label,trajectories in data:
    i = 0
    print "{:10s}: saving {:d} trajectories to readable format".format(label,len(trajectories))
    for traj in trajectories:
            filename = "{:s}-{:04d}.data".format(label,i)
            np.savetxt(filename,traj,fmt = '%.6e', delimiter = ' ')
            i += 1
