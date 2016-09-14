#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

def growth(n0,tau=6,a=1,y=2,s0=40000):
    return np.array([n*np.exp(a*tau) if n*np.exp(a*tau) < s0*y+n else s0*y+n for n in n0])


n0 = np.arange(1,100)

for dexp in range(1,5):
    dilution = np.power(10.,-dexp)
    print >> sys.stderr,dilution
    for sexp in np.arange(1,5,.5):
        print >> sys.stderr," "+str(sexp)
        n1 = dilution * growth(n0,s0=np.power(10.,sexp))
        
        for i in range(len(n0)):
            print dilution,sexp,n0[i],n1[i]

        print
