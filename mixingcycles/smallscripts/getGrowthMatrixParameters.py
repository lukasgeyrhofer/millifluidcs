#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc

def array_to_str(x):
    return ' '.join(['{:e}'.format(y) for y in x])


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
args = parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"could not open file"

if not g.hasGrowthMatrix():
    raise IOError,"pickle file does not contain a growth matrix"


out = ""

out += " -a " + array_to_str(g.growthrates)
out += " -y " + array_to_str(g.yieldfactors)

out += " -S {:e} ".format(g.env.substrate)

if not g.env.mixingtime is None:
    out += " -T {:e} ".format(g.env.mixingtime)

sys.stdout.write(out + '\n')

