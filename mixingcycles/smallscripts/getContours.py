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
parser.add_argument("-C","--contourvalues",default=[0],type=float,nargs="*")
parser.add_argument("-o","--outbasename",default=None)
args = parser.parse_args()


try:
    data = np.genfromtxt(args.infile)
except:
    raise IOError("could not open file '{}'".format(args.infile))

x = np.unique(data[:,args.xcolumn])
y = np.unique(data[:,args.ycolunm])
z = data[:,args.zcolumn].reshape((len(x),len(y)))

for cval in args.contourvalues:
    if not args.outbasename is None:
        fp = open(args.outbasename + "_{:.4e}.txt".format(cval),"w")
    else:
        fp = sys.stdout
        fp.write("# contourvalue = {:e}".format(cval))
        
    contours = measure.find_contours(z,cval)
    for c in contours:
        for i in range(len(c)):
            ix = int(np.floor(c[i,0]))
            iy = int(np.floor(c[i,1]))
            px = c[i,0] - ix
            py = c[i,1] - iy
            
            try:    cx = (1.-px)*x[ix] + px*x[ix+1]
            except: cx = x[ix]
                
            try:    cy = (1.-py)*y[iy] + py*y[iy+1]
            except: cy = y[iy]
            
            fp.write("{:.6f} {:.6f}\n".format(cx,cy))
        fp.write("\n")
    if not args.outbasename is None:
        fp.close()


