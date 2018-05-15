#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/..')

import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-x","--xcolumn",type=int,default=0)
parser.add_argument("-y","--ycolumn",type=int,default=1)
parser.add_argument("-z","--zcolumn",type=int,default=2)
parser.add_argument("-T","--transpose",action="store_true",default=False)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    raise IOError("could not open file '{}'".format(args.infile))

xlist = np.sort(np.unique(data[:,args.xcolumn]))
ylist = np.sort(np.unique(data[:,args.ycolumn]))
zdata = data[:,args.zcolumn].reshape((len(xlist),len(ylist)))

if args.transpose:
    zdata = zdata.T
    tmp   = xlist[:]
    xlist = ylist[:]
    ylist = tmp[:]

for i,x in enumerate(xlist):
    yi = zdata[i,:].argmax()
    ymax = ylist[yi]
    print "{} {} {}".format(x,ymax)
