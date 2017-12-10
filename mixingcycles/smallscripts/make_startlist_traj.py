#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest = "shape")

parser_rect = subparsers.add_parser("R")
parser_rect.add_argument("-x","--xrange",default=[10,80],type=float,nargs=2)
parser_rect.add_argument("-X","--xsteps",default=10,type=int)
parser_rect.add_argument("-y","--yrange",default=[10,60],type=float,nargs=2)
parser_rect.add_argument("-Y","--ysteps",default=10,type=int)

parser_circ = subparsers.add_parser("C")
parser_circ.add_argument("-r","--radius",default=10,type=float)
parser_circ.add_argument("-s","--steps",default=12,type=int)

args = parser.parse_args()

if args.shape == "R":
    x = np.linspace(start = args.xrange[0], stop = args.xrange[1], num = args.xsteps)
    y = np.linspace(start = args.yrange[0], stop = args.yrange[1], num = args.ysteps)

    nx = len(x)
    ny = len(y)

    startpos = np.concatenate([np.transpose([x,np.ones(nx)*y[0]]),np.transpose([np.ones(ny-2)*x[-1],y[1:-1]]),np.transpose([x[::-1],np.ones(nx)*y[-1]]),np.transpose([np.ones(ny-2)*x[0],y[1:-1]])])
elif args.shape == "C":
    assert args.radius > 0
    assert args.steps >= 2

    angles = 2 * math.pi * np.arange(args.steps)/(1. * args.steps)
    startpos = args.radius * np.transpose([np.sin(angles),np.cos(angles)])
else:
    raise NotImplementedError


for p in startpos:
    print '{} {}'.format(*p)



