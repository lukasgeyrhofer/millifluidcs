#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24)

args = parser.parse_args()


g = gc.GrowthDynamicsAntibiotics(**vars(args))









print g
