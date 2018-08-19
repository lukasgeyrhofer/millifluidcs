#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import itertools

def closeto(point1,point2,accepteddistance = 1e-3,restrictdimension = None):
    assert len(point1) == len(point2), "Coordinates must have identical dimension"
    if restrictdimension is None:
        return np.linalg.norm(point1-point2) < accepteddistance
    elif isinstance(restrictdimension,str):
        assert len(restrictdimension) == len(point1),"Restriction must have identical dimension as points"
        assert restrictdimension.count('1') + restrictdimension.count('0') == len(restrictdimension),"Restriction on dimensions has to be either string of '1' and '0' or a single number"
        d = 0.
        for i in range(len(restrictdimension)):
            if restrictdimension[i] == '1':
                d += (point1[i] - point2[i])**2
        return np.sqrt(d) < accepteddistance
    elif isinstance(restrictdimension,int):
        assert 0 <= restrictdimension < len(point1),"Restriction must chose valid dimension"
        return np.sqrt(point1[restrictdimension] - point2[restrictdimension]) < accepteddistance


def invert(bitstring):
    return "".join([str(1-int(x)) for x in bitstring])


def writeoutput(filename,trajectories):
    fp = open(filename,"w")
    for t in trajectories:
        for p in t:
            print >> fp,"".join([" {}".format(x) for x in p])
        print >> fp
    fp.close()

parser = argparse.ArgumentParser()
parser.add_argument("-i","--trajectoryfile")
parser.add_argument("-d","--accepteddistance",type=float,default=1e-3)
parser.add_argument("-o","--basenameoutfile",default=None)
args = parser.parse_args()

trajectories = list()
current      = list()
dim          = 0
try:
    fp = open(args.trajectoryfile,"r")
except:
    raise IOError,"Error opening file '{}'".format(args.trajectoryfile)


for line in fp.readlines():
    v = line.split()
    if len(v) < 2:
        # usually an empty line ...
        trajectories.append(np.array(current))
        current = list()
    else:
        p = np.array(v,dtype=np.float)
        dim = len(p)
        current.append(p)
fp.close()


# create strings of 0's and 1's, where 0 means trajectory reached zero in population (indicated by position of 0 in the string)
zero                       = np.zeros(dim)
fixedpoint_restrictions    = sorted([''.join(x) for x in itertools.product('01',repeat = dim)],cmp = lambda a,b: np.sign(b.count('0') - a.count('0')))
sortedtrajectories         = dict()
for fp in fixedpoint_restrictions:
    sortedtrajectories[fp] = list()

for t in trajectories:
    for fp in fixedpoint_restrictions:
        # index 1 in the 'closeto' function means use and compare this dimension
        # index 0 in the output means that the endpoint of the trajectory is zero in this dimension (while 1 indicates endpoint is away from boundary in this dimension)
        # both contradict each other, thus invert the bitstring!
        if closeto(t[-1],zero,restrictdimension = invert(fp)):
            sortedtrajectories[fp].append(t)
            break

if args.basenameoutfile is None:    outfile = args.trajectoryfile
else:                               outfile = args.basenameoutfile

for fp in fixedpoint_restrictions:
    writeoutput(outfile + "_" + fp,sortedtrajectories[fp])

