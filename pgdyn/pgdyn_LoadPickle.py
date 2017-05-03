#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-i","--picklefile")
args = parser.parse_args()


g = pickle.load(open(args.picklefile))

g.ExtendGrowthMatrix(5)


print g.growthmatrix
