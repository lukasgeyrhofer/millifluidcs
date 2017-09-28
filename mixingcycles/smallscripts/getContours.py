#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import matplotlib.pyplot as plt
from skimage import measure

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-x","--xcolumn",default=0,type=int)
parser.add_argument("-y","--ycolunm",default=1,type=int)
parser.add_argument("-z","--zcolumn",default=2,type=int)
parser.add_argument("-C","--contourvalue",default=0,type=float)
args = parser.parse_args()


try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

x = np.unique(data[:,args.xcolumn])
y = np.unique(data[:,args.ycolunm])
z = data[:,args.zcolumn].reshape((len(x),len(y)))
contours = measure.find_contours(z,args.contourvalue)

for c in contours:
    for i in range(len(c)):
        ix = int(np.floor(c[i,0]))
        iy = int(np.floor(c[i,1]))
        px = c[i,0] - ix
        py = c[i,1] - iy
        
        try:
            cx = (1.-px)*x[ix] + px*x[ix+1]
        except:
            cx = x[ix]
            
        try:
            cy = (1.-py)*y[iy] + py*y[iy+1]
        except:
            cy = y[iy]
        
        print "{:.6f} {:.6f}".format(cx,cy)
    print


