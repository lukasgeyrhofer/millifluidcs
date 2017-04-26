#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-T","--threshold",type=float,default=3e-2)
args = parser.parse_args()

data = dict()

for filename in args.infiles:
    print filename
    try:
        data[filename] = np.genfromtxt(filename,delimiter=',',names = True, dtype=float)
        newfilename = filename.replace("csv","data")

        t = data[filenale]['time']
        b = data[filename]['Channel1_mean']

        
        #print t,b
        fp = open(newfilename,"w")
        for dt in t:
            print >> fp, dt
        fp.close()
    except:
        continue
