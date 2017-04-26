#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
args = parser.parse_args()

for filename in args.infiles:
    try:
        data = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
        newdata = np.transpose(np.array([data['time'],data['Channel1_mean']]))
        np.savetxt(filename.replace("csv","data"),newdata)
    except:
        continue
