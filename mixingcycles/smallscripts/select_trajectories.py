#!/usr/bin/env python

import argparse
import numpy as np
import math,sys

def distance(x,y):
    return np.sqrt(np.dot(x-y,x-y))

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",required=True)
parser.add_argument("-n","--coordinates",nargs=2,type=float,default=[1,1])
parser.add_argument("-r","--radius",default=.1,type=float)
parser_mode = parser.add_mutually_exclusive_group(required=True)
parser_mode.add_argument("-S","--start",action="store_true",default=False)
parser_mode.add_argument("-P","--passthrough",action="store_true",default=False)
args = parser.parse_args()

try:
    fp = open(args.infile)
except:
    raise IOError("could not open file '{}'".format(args.infile))

coordinates  = np.array(args.coordinates,dtype=np.float)
trajectories = list()

start_reading = True
accept_traj   = 0
cur_traj      = list()

for line in fp.readlines():
    values = np.array(line.split(),dtype=np.float)
    if len(values) == 0:
        if accept_traj > 0:
            trajectories.append(np.array(cur_traj))

        start_reading = True
        cur_traj      = list()
        accept_traj   = 0
    else:
        if args.passthrough or (args.start and start_reading):
            if distance(values,coordinates) < args.radius:
                accept_traj += 1

        cur_traj.append(values)
        start_reading = False
        
for t in trajectories:
    for pos in t:
        print '{:.6e} {:.6e}'.format(pos[0],pos[1])
    print 



