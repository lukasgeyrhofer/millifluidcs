#!/usr/bin/env python


import numpy as np
import argparse
import sys

sys.path.append(sys.path[0] + '/../mixingcycles')
import growthclasses as gc


parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)

parser_PG = parser.add_argument_group(description = "Parameters for interactions with public good")
parser_PG.add_argument("-P","--pgproduction",nargs="*",default=np.zeros(2))
parser_PG.add_argument("-A","--pginteractiongrowthrates",nargs="*",default=np.zeros(2))
parser_PG.add_argument("-Y","--pginteractionyieldfactor",nargs="*",default=np.zeros(2))
args = parser.parse_args()

g = gc.GrowthDynamicsPublicGoods(**vars(args))


