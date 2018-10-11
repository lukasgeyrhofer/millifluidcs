#!/usr/bin/env python3
#-*- coding: utf-8 -*-

'''
==================================
=  mixingcycles_SingleStep.py
==================================

    Computes a single step of average inoculum sizes over multiple cycles


    Lukas Geyrhofer, l.geyrhofer@technion.ac.il, 2018

'''


import argparse
import numpy as np
import sys,math
import pickle
import growthclasses as gc
import itertools


def main():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infile",required = True)
    parser_io.add_argument("-o","--outfile",required = True)
    parser_io.add_argument("-v","--verbose",default=False,action="store_true")

    parser = gc.AddDilutionParameters(parser)
    parser = gc.AddLatticeParameters(parser)
    
    args = parser.parse_args()

    g           = gc.LoadGM(**vars(args))
    dlist       = gc.getDilutionList(**vars(args))
    axis1,axis2 = gc.getInoculumAxes(**vars(args)) # either (n,x) or [ (n1,n2) if args.AbsoluteCoordinates == True ]

    gm1 = g.growthmatrix[:,:,0]
    gm2 = g.growthmatrix[:,:,1]

    for dilution in dlist:
        if args.verbose:
            sys.stderr.write("# computing single step dynamics for D = {:e}\n".format(dilution))
            
        fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")
        for i,a1 in enumerate(axis1):
            for j,a2 in enumerate(axis2):
                next1 = gc.SeedingAverage(gm1, gc.TransformInoculum([a1,a2],args.AbsoluteCoordinates,True)) * dilution
                next2 = gc.SeedingAverage(gm2, gc.TransformInoculum([a1,a2],args.AbsoluteCoordinates,True)) * dilution
                fp.write('{} {} {} {}\n'.format(a1,a2,*gc.TransformInoculum([next1,next2],True,args.AbsoluteCoordinates)))
            fp.write('\n')
        fp.close()
                
                
if __name__ == "__main__":
    main()
