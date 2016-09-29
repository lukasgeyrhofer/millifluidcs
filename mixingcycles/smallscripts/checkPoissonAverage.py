#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + "/..")
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

mlogm = m
mlogm[m>0] *= np.log(m[m>0])

for aval in a:
    for nval in n:
        px = gc.PoissonSeedingVectors(m,[nval])
        ma = m**aval
        ma[ma is np.nan] = 0
        
        expo = np.dot(ma,px[0])
        eps1 = np.dot(m + (aval-1)*mlogm,px[0])
        eps0 = np.dot(1 + (aval)*np.log(m[1:]),px[0][1:])
        
        apprexpo = nval**aval
        appreps1 = nval*(1+(aval-1)*np.log(nval))
        appreps0 = 1 + (aval)*np.log(nval)
        
        print(aval,nval,expo,eps1,eps0,apprexpo,appreps1,appreps0)
    print()



