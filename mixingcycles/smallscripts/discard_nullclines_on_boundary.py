#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile",default=None)
args = parser.parse_args()

try:
    fp = open(args.infile,"r")
except:
    raise IOError("could not open file '{}'".format(args.infile))


if args.outfile is None:
    out = sys.stdout
else:
    out = open(args.outfile,"w")

cur_traj = list()
print_traj = True
for line in fp.readlines():
    values = np.array(line.split(),dtype=np.float)
    if len(values) == 0 and len(curtraj) > 0:
        if print_traj:
            for a,b in cur_traj:
                out.write("{} {}\n".format(a,b))
            out.write("\n")
        cur_traj = list()
        print_traj = True
    else:
        cur_traj.append(values)
        if values[0] == 0 or values[1] == 0:
            print_traj = False

if not args.outfile is None:
    out.close()


