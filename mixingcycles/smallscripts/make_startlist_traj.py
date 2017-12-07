#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-x","--xrange",default=[10,80],type=float,nargs=2)
parser.add_argument("-X","--xstep",default=10,type=float)
parser.add_argument("-y","--yrange",default=[10,60],type=float,nargs=2)
parser.add_argument("-Y","--ystep",default=10,type=float)
args = parser.parse_args()

x = np.arange(start = args.xrange[0], stop = args.xrange[1], step = args.xstep)
y = np.arange(start = args.yrange[0], stop = args.yrange[1], step = args.ystep)

nx = len(x)
ny = len(y)

startpos = np.concatenate([np.transpose([x,np.ones(nx)*y[0]]),np.transpose([np.ones(ny-2)*x[-1],y[1:-1]]),np.transpose([x,np.ones(nx)*y[-1]]),np.transpose([np.ones(ny-2)*x[0],y[1:-1]])])

for p in startpos:
    print '{} {}'.format(*p)



