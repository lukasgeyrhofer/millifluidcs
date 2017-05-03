#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
ioparser.add_argument("-i","--infiles",nargs="*")
ioparser.add_argument("-o","--outbasename",default = None)
ioparser.add_argument("-t","--templatefile",default = None)
ioparser.add_argument("-r","--restrictionfile",default = None)

dparser = parser.add_argument_group(description = "==== Data parameters ====")
dparser.add_argument("-C","--columns",nargs="*", default = ["Channel1_mean"])
dparser.add_argument("-c","--signalcolumn",default = 0,type = float)
dparser.add_argument("-B","--splitBackForthTrajectories",default = False,action="store_true")
dparser.add_argument("-x","--lowertimecutoff",default = None,type = float)
dparser.add_argument("-X","--uppertimecutoff",default = None,type = float)
dparser.add_argument("-s","--lowersignalcutoff",default = None,type = float)
dparser.add_argument("-S","--uppersignalcutoff",default = None,type = float)
args = parser.parse_args()

if not "time" in args.columns:
    columns = ["time"] + args.columns
    timecolumnoffset = 1
else:
    columns = args.columns
    timecolumnoffset = 0



data = mdc.DropletData(infiles = args.infiles, templatefile = args.templatefile, datacolumns = columns, splitBackForthTrajectories = True)
if not args.lowertimecutoff is None:
    data.set_restriction("time","min",args.lowertimecutoff)
if not args.uppertimecutoff is None:
    data.set_restriction("time","max",args.uppertimecutoff)
if not args.lowersignalcutoff is None:
    data.set_restriction(columns[args.signalcolumn + timecolumnoffset],"min",args.lowersignalcutoff)
if not args.uppersignalcutoff is None:
    data.set_restriction(columns[args.signalcolumn + timecolumnoffset],"max",args.uppersignalcutoff)

if not args.restrictionfile is None:
    data.load_restrictions_from_file(args.restrictionfile)

for label,trajectories in data:
    i = 0
    for traj in trajectories:
        if len(traj) > 0:
            if not args.outbasename is None:
                filename = args.outbasename + "{:s}-{:04d}.data".format(label,i)
            else:
                filename = "{:s}-{:04d}.data".format(label,i)
            np.savetxt(filename,traj,fmt = '%.6e', delimiter = ' ')
            i += 1
    print "{:12s}: saving {:4d} trajectories to readable format".format(label,i)
