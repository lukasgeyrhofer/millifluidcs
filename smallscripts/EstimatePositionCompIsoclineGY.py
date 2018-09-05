#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

from scipy.stats import poisson


def func(n,p,s):
    return -p/n + (1.-1./n)*np.log(s/n)


parser = argparse.ArgumentParser()
parser.add_argument("-p","--parameterratio",default=1,type=float)
parser.add_argument("-S","--resources",default=1e5,type=float)
parser.add_argument("-m","--maxn",default=100,type=int)

parser.add_argument("-M","--maxSteps",default=1000,type=int)
parser.add_argument("-P","--precision2",default=1e-20,type=float)
parser.add_argument("-A","--alpha",default=1.,type=float)
parser.add_argument("-q","--quiet",default=False,action="store_true")
args=  parser.parse_args()

m    = np.arange(start = 1,stop = args.maxn)
n0   = 0
n1   = args.parameterratio/np.log(args.resources)
step = 0


while True:
    if ((n1-n0)**2 < n1**2*args.precision2) or (step > args.maxSteps):
        break
    
    n0    = n1
    
    p     = poisson.pmf(m,n0)
    f     = np.dot(p,func(m,args.parameterratio,args.resources))
    fp    = np.dot(p,func(m+1,args.parameterratio,args.resources) - func(m,args.parameterratio,args.resources))
    
    n1   -= args.alpha * f/fp
    step += 1

if args.quiet:
    sys.stdout.write("{:.6e}\n".format(n1))
else:
    sys.stdout.write("{:6e} {:6e} {:.6e} {:d}\n".format(args.parameterratio,args.resources,n1,step))
