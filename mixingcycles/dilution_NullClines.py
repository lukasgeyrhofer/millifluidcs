#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import matplotlib.pyplot as plt
from skimage import measure

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",nargs="*")
args = parser.parse_args()


for fn in args.infile:
    print fn
    try:
        data = np.genfromtxt(fn)
    except:
        print >> sys.stderr,"could not open file"
        exit(1)


    x = np.unique(data[:,0])
    y = np.unique(data[:,1])
    n1 = data[:,2].reshape((len(x),len(y)))
    n2 = data[:,3].reshape((len(x),len(y)))

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    contours1 = measure.find_contours(n1-np.outer(x,np.ones(len(y))),0.)
    contours2 = measure.find_contours(n2-np.outer(np.ones(len(x)),y),0.)


    fp1 = open(fn+".n1","w")
    for c in contours1:
        for i in range(len(c)):
            print >> fp1,c[i,0]*dx,c[i,1]*dy
        print >> fp1
    fp1.close()

    fp2 = open(fn+".n2","w")
    for c in contours2:
        for i in range(len(c)):
            print >> fp2,c[i,0]*dx,c[i,1]*dy
        print >> fp2
    fp2.close()

