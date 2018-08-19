#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

def sigmoid(x,m,s):
    return 1/(1+np.exp(-(x-m)/s))



parser = argparse.ArgumentParser()
parser.add_argument("-k","--numdroplets",type=int,default=1000)
parser.add_argument("-p","--periodsecond",type=int,default=200)
args = parser.parse_args()



k = np.arange(args.numdroplets)

x1 = sigmoid(k,0.5*args.numdroplets,0.1*args.numdroplets)
x2 = np.zeros(args.numdroplets)
start = 0
while start <= args.numdroplets:
    x2[start:start + args.periodsecond/2] = sigmoid(k[start:start + args.periodsecond/2],start + 0.25*args.periodsecond,0.05*args.periodsecond)
    x2[start + args.periodsecond/2:start + args.periodsecond] = 1 - sigmoid(k[start + args.periodsecond/2:start + args.periodsecond],start + 0.75*args.periodsecond,0.05*args.periodsecond)
    start += args.periodsecond

for a in zip(k,x1,x2):
    print "{} {} {}".format(*a)







