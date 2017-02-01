#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


def dist(p1,listp2):
    return np.array([np.linalg.norm(np.log(p1)-np.log(p2)) for p2 in listp2])


parser = argparse.ArgumentParser()
parser.add_argument("-1","--infile1")
parser.add_argument("-2","--infile2")
parser.add_argument("-i","--interpolate",type=float,default=.5)
args = parser.parse_args()


try:
    data1 = np.genfromtxt(args.infile1)
    data2 = np.genfromtxt(args.infile2)
except:
    print >> sys.stderr,"could not open file"
    exit(1)


for p in data1:
    newp = args.interpolate * p + (1-args.interpolate) * data2[dist(p,data2).argmin()]
    print "{:.6} {:.6}".format(*newp)
