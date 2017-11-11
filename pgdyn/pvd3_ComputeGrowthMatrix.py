#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/../mixingcycles/')
import growthclasses as gc

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser,dilution=True)

parser_pvd = parser.add_argument_group(description = "==== Parameters for interactions with PVD ====")
parser_pvd.add_argument("-Y","--PVD_Internal_Yield",type=float,nargs="*",default=[1,1])
parser_pvd.add_argument("-P","--PVD_Production_Efficiency",type=float,nargs="*",default=[1e-3,0])
parser_pvd.add_argument("-M","--PVD_Matching_Receptors",type=float,nargs="*",default=[1,1])
parser_pvd.add_argument("-I","--PVD_Initial_Internal_Iron",type=float,nargs="*",default=[0,0])

