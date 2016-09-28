#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math


import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-M","--maxM",type=int,default=300)
parser.add_argument("-N","--maxN",type=float,default=10)
parser.add_argument("-n","--dN",type=float,default=.1)
parser.add_argument("-A","--maxA",type=float,default=2)
parser.add_argument("-a","--dA",type=float,default=.1)
args = parser.parse_args()


a = np.arange(-args.maxA,args.maxA+args.dA,args.dA)
n = np.arange(int(np.ceil(args.maxN/args.dN))+1)*args.dN
m = np.arange(args.maxM)


for aval in a:
    for nval in n:
        px = gc.PoissonSeedingVectors(m,[nval])
        ma = m**aval
        ma[ma is np.nan] = 0
        print(aval,nval,np.dot(ma,px[0]),nval**aval)
    print()



