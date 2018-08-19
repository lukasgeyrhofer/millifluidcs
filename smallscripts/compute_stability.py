#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

def metric(x):
    return np.sqrt(np.dot(x,x))

def getStability(x):
    unstable = np.sum(np.where([metric(x[:2]) >= 1, metric(x[2:4]) >= 1 ],1,0))
    if metric(x[1]) > 0:
        unstable += 3
    return unstable


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",default=None)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-F","--firstcolumn",default=None,type=float)
args = parser.parse_args()

if not args.infile is None:
    try:
        data = np.genfromtxt(args.infile)
    except:
        raise IOError
else:
    raw = list()
    for line in sys.stdin.readline():
        raw.append(line)
    try:
        data = np.array(raw,dtype=np.float)
    except:
        raise ValueError
    
out = np.zeros(len(data))

for i,val in enumerate(data):
    if (val[3] > 0) and (val[4] > 0):
        ev = np.array([val[6],val[8],val[7],val[9]],dtype=np.float)
        out[i] = getStability(ev)
    else:
        out[i] = -1

if not args.firstcolumn is None:
    fc = np.ones(len(out)) * args.firstcolumn
    outdata = np.transpose([fc,data[:,0],out])
else:
    outdata = np.transpose([data[:,0],out])

if not args.outfile is None:
    np.savetxt(args.outfile,outdata)
else:
    for x in outdata:
        for val in x:
            sys.stdout.write(" {:6e}".format(val))
        sys.stdout.write("\n")
