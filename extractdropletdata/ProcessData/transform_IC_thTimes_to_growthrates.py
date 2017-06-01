#!/usr/bin/env python

import argparse
import numpy as np
import sys,math


def growthrate(time,n0,n1):
    return np.log(n1/n0)/time


parser = argparse.ArgumentParser()
parser.add_argument("-n","--initialconditionsPoisson",default=25,type=int)
parser.add_argument("-N","--Nmax",default=100,type=int)
parser.add_argument("-t","--thresholdtimefile",default=None)
parser.add_argument("-T","--threshold",type=float,default=200)

parser.add_argument("-a","--amin",default=0,type=float)
parser.add_argument("-A","--amax",default=2,type=float)
parser.add_argument("-d","--da",default=.1,type=float)
args = parser.parse_args()


n  = np.arange(args.Nmax)
pn = np.exp(-args.initialconditionsPoisson)*np.power(args.initialconditionsPoisson,n)/np.gamma(n-1)

a  = np.arange(start = args.amin,stop = args.amax, step = args.da)
pa = np.zeros(len(a))

try:
    thdata = np.genfromtxt(args.thresholdtimefile)
    t  = thdata[:,0]
    pt = thdata[:,1]
except:
    raise IOError

print pn





