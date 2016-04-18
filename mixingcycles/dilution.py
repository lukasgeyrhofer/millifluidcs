#!/usr/bin/env python

import numpy
import argparse
import sys,math
from scipy.misc import factorial

def get_growth(initialcells,**parameters):
    r = np.zeros(len(initialcells))
    
    
    
    
    
    return r



parser = argparse.ArgumentParser()
parser.add_argument("-n","--maxsize",type=int,default=100)
args = parser.parse_args()

p = {}
p['yield'] = np.array([1,2])
p['growth'] = np.array([2,1])
p['substrate'] = 1e4
p['mixingtime'] = 6
p['dilution'] = 2e-4


n0 = np.zeros((args.maxsize,args.maxsize,2))
for i in range(args.maxsize):
    for j in range(args.maxsize):
        n0[i,j] = get_growth(np.array([i,j]),p)


n1 = np.zeros((args.maxsize,args.maxsize,2))
x = np.arange(args.maxsize)
fx = factorial(x,exact=False)

px = np.zeros((args.maxsize,args.maxsize))
for i in range(args.maxsize):
    px[i] = np.power(i,x)*np.exp(-i)/fx

for i in range(args.maxsize):
    for j in range(args.maxsize):
        n1[i,j] = np.dot(np.dot(n0,px[i]),px[j])
        
        print i,j,n1[i,j,0],n1[i,j,1]
    print













