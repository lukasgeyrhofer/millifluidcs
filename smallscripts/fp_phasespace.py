#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
args = parser.parse_args()


try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

a_all  = data[:,0]
y_all  = data[:,1]
n1_all = data[:,2]
n2_all = data[:,3]


alist = np.unique(a_all)
lastc = 0
for a in alist:
    y  = y_all [a_all == a]
    c = len(y)

    if c >= lastc:
        
        n1 = n1_all[a_all == a]
        n2 = n2_all[a_all == a]

        #print n1,n2

        n1r = n1[n1 > 0]
        y1r = y [n1 > 0]
        
        n2r = n2[n2 > 0]
        y2r = y [n2 > 0]
        
        try:
            # trivial approximation
            # y1 = y1r[-1]
            
            # linear interpolation
            dy1 = y1r[-1] - y1r[-2]
            dn1 = n1r[-1] - n1r[-2]
            y1  = y1r[-1] + n1r[-1]/dn1 * dy1
        except:
            # lists are not long enough, etc ...
            y1 = np.nan
        try:
            # trivial approximation
            # y2 = y2r[0]
            
            # linear interpolation
            dn2 = n2r[0] - n2r[1]
            dy2 = y2r[0] - y2r[1]
            y2 = y2r[0] - n2r[0]/dn2 * dy2
        except:
            y2 = np.nan
        
        
        print "{:.6f} {:.6f} {:.6f}".format(a,y1,y2)

    lastc = c






