#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math
sys.path.append(sys.path[0] + '/..')

import growthclasses as gc
from scipy.stats import poisson
import itertools

def func(ax1,ax2,params):
    x = itertools.product(ax1,ax2)
    print(x)
    return None


def main():

    parser = argparse.ArgumentParser()
    parser_param = parser.add_argument_group(description = "==== Parameters ====")
    parser_param.add_argument("-s","--yieldinc",default=2,type=float)
    parser_param.add_argument("-c","--expconst",default=1,type=float)
    parser_param.add_argument("-m","--maxpoisson",default=100,type=int)

    parser = gc.AddLatticeParameters(parser)

    args = parser.parse_args()

    axis1,axis2 = gc.getInoculumAxes(**vars(args)) # either (n,x) or [ (n1,n2) if args.AbsoluteCoordinates == True ]
    shape       = (len(axis1),len(axis2))
    m           = np.arange(args.maxpoisson)
    xw1 = np.zeros(shape)

    params = {    'sigma':      args.yieldinc,
                  'expconst':   args.expconst }

    for i,a1 in enumerate(axis1):
        for j,a2 in enumerate(axis2):
            inoc = gc.TransformInoculum([a1,a2],inabs = args.AbsoluteCoordinates, outabs = True)
            xw1[i,j] = gc.SeedingAverage(func(axis1,axis2,params), inoc)
            
            print("{} {} {}".format(a1,a2,xw1[i,j]))
        print("")

if __name__ == "__main__":
    main()






