#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-r","--radius",default=10,type=float)
parser.add_argument("-s","--steps",default=12,type=int)
args=parser.parse_args()

assert args.radius > 0
assert args.steps >= 2

angles = 2 * math.pi * np.arange(args.steps)/(1. * args.steps)


x = args.radius * np.transpose([np.sin(angles),np.cos(angles)])

for y in x:
    print "{:13.6e} {:13.6e}".format(*y)
