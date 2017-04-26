#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-L","--forcelen",type=int,default=None)
args = parser.parse_args()

for filename in args.infiles:
    try:
        data = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
        newdata = np.transpose(np.array([data['time'],data['Channel1_mean']]))
        if not args.forcelen is None:
            if len(newdata) < args.forcelen:
                zeros = np.zeros(shape = (args.forcelen- len(newdata),2))
                newdata = np.concatenate((newdata,zeros))
                print "forcelength:",filename
        np.savetxt(filename.replace("csv","data"),newdata,fmt='%.6e')
    except:
        continue
