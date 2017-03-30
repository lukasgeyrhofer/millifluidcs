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
        assert restrictdimension.count('1') > 0,"Restriction must contain at least one dimension to compare"
        d = 0.
        for i in range(len(restrictdimension)):
            if restrictdimension[i] == '1':
                d += (point1[i] - point2[i])**2
        return np.sqrt(d) < accepteddistance
    elif isinstance(restrictdimension,int):
        assert 0 <= restrictdimension < len(point1),"Restriction must chose valid dimension"
        return np.sqrt(point1[restrictdimension] - point2[restrictdimension]) < accepteddistance


def writeoutput(filename,trajectories):
    fp = open(filename,"w")
    for t in trajectories:
        for p in t:
            print >> fp,p[0],p[1]
        print >> fp
    fp.close()

parser = argparse.ArgumentParser()
parser.add_argument("-i","--trajectoryfile")
parser.add_argument("-d","--accepteddistance",type=float,default=1e-3)
args = parser.parse_args()

trajectories = list()
current      = list()
dim          = 0
try:
    fp = open(args.trajectoryfile)
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
fixedpoint_restrictions = sorted([''.join(x) for x in itertools.product('01',repeat = dim)],cmp = lambda a,b: np.sign(b.count('0') - a.count('0')))
sortedtrajectories      = dict()

for t in trajectories:
    for fp in fixedpoint_restrictions:
        if not sortedtrajectories.has_key(fp):
            sortedtrajectories[fp] = list()
        if closeto(t[-1],zero,fp):
            sortedtrajectories.append(t)
            break

for fp in fixedpoint_restrictions:
    writeoutput(args.infile + "_" + fp,sortedtrajectories[fp])

