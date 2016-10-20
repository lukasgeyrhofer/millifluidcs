#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

from skimage import measure


def write_contours(contours,filename,y,a):
    fp = open(filename,"w")
    for c in contours:
        for i in range(len(c)):
            print >> fp,y[c[i,0]],a[c[i,1]]
        print >> fp
    fp.close()


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-d","--dp",type=float,default=.1)
args = parser.parse_args()



try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

yvalues = np.unique(data[:,0])
avalues = np.unique(data[:,1])

y =              np.repeat([yvalues],len(avalues),axis=0)
a = np.transpose(np.repeat([avalues],len(yvalues),axis=0))



gamma = (data[:,2]).reshape(np.shape(y))
G = 1 - 2*gamma


p = args.dp
while p <= 1:
    print p
    stab0 = 1./(y*(1.-p*G))
    stab1 = y/(1.+p*G)
    stabi = 2+(1-y)*(2-(1+y)*(1-p*G))/(y*p*p*G*G)
    fpi   = y/(y-1.) - 1/(p*G)

    try:
        contours0p = measure.find_contours(stab0, 1.)
        contours0m = measure.find_contours(stab0,-1.)
        write_contours(contours0p,"cont.stab0.+1.p%5.3f.txt"%p,yvalues,avalues)
        write_contours(contours0m,"cont.stab0.-1.p%5.3f.txt"%p,yvalues,avalues)
    except:
        pass

    try:
        contours1p = measure.find_contours(stab1, 1.)
        contours1m = measure.find_contours(stab1,-1.)
        write_contours(contours1p,"cont.stab1.+1.p%5.3f.txt"%p,yvalues,avalues)
        write_contours(contours1m,"cont.stab1.-1.p%5.3f.txt"%p,yvalues,avalues)
    except:
        pass

    try:
        contoursip = measure.find_contours(stabi, 1.)
        contoursim = measure.find_contours(stabi,-1.)
        write_contours(contoursip,"cont.stabi.+1.p%5.3f.txt"%p,yvalues,avalues)
        write_contours(contoursim,"cont.stabi.-1.p%5.3f.txt"%p,yvalues,avalues)
    except:
        pass

    try:
        contoursFP0 = measure.find_contours(fpi,0.)
        contoursFP1 = measure.find_contours(fpi,1.)
        write_contours(contoursFP0,"cont.fpi.0.p%5.3f.txt"%p,yvalues,avalues)
        write_contours(contoursFP1,"cont.fpi.1.p%5.3f.txt"%p,yvalues,avalues)
    except:
        pass

    p += args.dp
