#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

def mymin(x):
    if len(x) > 0:
        return np.min(x)
    else:
        return None
def mymax(x):
    if len(x) > 0:
        return np.max(x)
    else:
        return None
    


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
args = parser.parse_args()


try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

a = data[:,0]
y = data[:,1]
fp1 = data[:,2]
fp2 = data[:,3]


alist = np.unique(a)

for aa in alist:
    yr   = y[a == aa]
    fpr1 = fp1[a == aa]
    fpr2 = fp2[a == aa]

    y1 = mymin(yr[fpr1 == 0])
    y2 = mymax(yr[fpr2 == 0])
    
    print aa,y1,y2








