#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,defaultmixingtime=24,dilution=True)

parser_ab = parser.add_argument_group(description = "Parameters for dynamics of PVD")
parser_ab.add_argument("-V","--PVDconc",type=float,default=0.0)
parser_ab.add_argument("-k","--PVDincreaseS",type=float,default=.2)
parser_ab.add_argument("-K","--PVDmaxFactorS",type=float,default=1.2)
parser_ab.add_argument("-p","--PVDproduction",nargs="*",default=[1,0])

parser_iterationmap = parser.add_argument_group(description = "Parameters for iterationmap")
parser_iterationmap.add_argument("-m","--maxsize",type=int,default=100)
parser_iterationmap.add_argument("-M","--step",type=int,default=1)
parser_iterationmap.add_argument("-n","--outputmax",type=float,default=20)
parser_iterationmap.add_argument("-d","--outputdx",type=float,default=.1)
parser_iterationmap.add_argument("-P","--poissonseeding",default=False,action="store_true")

parser.add_argument("-o","--outfile",default=None)

args = parser.parse_args()
g = gc.GrowthDynamicsPyoverdin(**vars(args))


print g.ParameterString()
