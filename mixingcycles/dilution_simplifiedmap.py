#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("-n","--etamin",type=float,default=1)
parser.add_argument("-N","--etamax",type=float,default=50)
parser.add_argument("-d","--etastep",type=float,default=.1)
parser.add_argument("-a","--alphamin",type=float,default=1e-2)
parser.add_argument("-A","--alphamax",type=float,default=1e2)
parser.add_argument("-D","--alphastep",type=float,default=.1)

pa = parser.add_argument_group(description="algorithm parameters")
pa.add_argument("-m","--maxsteps",type=int,default=10000)
pa.add_argument("-p","--precision",type=float,default=1e-10)
pa.add_argument("-R","--NRalpha",type=float,default=1)

args = parser.parse_args()

alist   = np.exp(np.arange(np.log(args.alphamin),np.log(args.alphamax),args.alphastep))
etalist = np.exp(np.arange(np.log(args.etamin),np.log(args.etamax),args.etastep))

prec2 = args.precision**2

f = np.zeros(2)
j = np.zeros((2,2))

for a,e in itertools.product(alist,etalist):
    n = np.ones(2) * (a + np.log(e))/2. # mean of n ~ O(a)
                                        # and     n ~ O(exp(eta))
                                        # as starting value
    dn = n
    i = 0
    # NR iterations to compute fixed point numerically
    while np.sum((dn/n)**2) > prec2:        
        #print i,n
        # first evaluation of function f(n)
        f[0] = a * (1-np.exp(-n[0])) - n[0]
        f[1] = a * np.exp(-n[0]) * (1-np.exp(-n[1])) * e - n[1]
        
        # then computation of Jacobian J at the current estimate of the fixed point
        j[0,0] = a * np.exp(-n[0]) - 1
        j[0,1] = 0
        j[1,0] = - f[1]
        j[1,1] = a * np.exp(-n[0]-n[1]) * e - 1
        
        # invert Jacobian
        iJ = np.linalg.inv(j)
        
        
        # NR step
        dn  = -args.NRalpha * np.dot(iJ,f)
        n  += dn
        
        #print j,iJ,f,n,dn
        # number of cells in inoculum cannot be negative
        n[n<0]=0
        
        i += 1
        
        if i > args.maxsteps:
            break
        
    print "{:.6e} {:.6e} {:.6e} {:.6e} {:6d}".format(a,e,n[0],n[1],i)





















