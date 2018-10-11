#!/usr/bin/env python3
#-*- coding: utf-8 -*-

'''
==================================
=  mixingcycles_Trajectories.py
==================================

    Computes trajectories of average inoculum sizes over multiple cycles


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
    
    parser_lattice = parser.add_argument_group(description = "==== Lattice parameters ====")
    parser_lattice.add_argument("-A","--AbsoluteCoordinates",default=False,action="store_true",help="Use (n1,n2) instead of (n,x) as coordinates")
    parser_lattice_startingconditions = parser_lattice.add_mutually_exclusive_group()
    parser_lattice_startingconditions.add_argument("-N","--maxInoculum",type=float,default=40)
    parser_lattice_startingconditions.add_argument("-I","--initialcoordinatesfile",default=None)
    parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
    parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)

    args = parser.parse_args()


    g     = gc.LoadGM(**vars(args))
    dlist = gc.getDilutionList(**vars(args))

    if args.initialcoordinatesfile is None:
        axis1,axis2 = gc.getInoculumAxes(**vars(kwargs))
        coordinates = list(itertools.product(axis1,axis2))
    else:
        try:
            fp_coords = open(args.initialcoordinatesfile)
        except:
            raise IOError("could not open file '{}' to load coordinates".format(args.initialcoordinatesfile))
        initialcoordinates = list()
        for line in fp_coords.readlines():
            try:
                values = np.array(line.split(),dtype=np.float)
                if len(values) >= 2:
                    initialcoordinates.append(gc.getAbsoluteInoculumNumbers(values[:2],args.AbsoluteCoordinates))
            except:
                continue
            fp_coords.close()


    mx,my = g.growthmatrixgrid
    gm1   = g.growthmatrix[:,:,0]
    gm2   = g.growthmatrix[:,:,1]

    if args.verbose:
        sys.stderr.write(str(g))

    for dilution in dlist:
        if args.verbose:
            sys.stderr.write("# computing trajectories for D = {:e}\n".format(dilution))
        fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")
        for ic1,ic2 in initialcoordinates:
            n1 = ic1
            n2 = ic2
            fp.write("{} {}\n".format(*gc.getCoordinatesFromAbsoluteInoculum([n1,n2],args.AbsoluteCoordinates)))
            for i in range(args.trajectorylength):
                next1 = gc.SeedingAverage(gm1,[n1,n2]) * dilution
                next2 = gc.SeedingAverage(gm2,[n1,n2]) * dilution
                n1,n2 = next1,next2
                fp.write("{} {}\n".format(*gc.getCoordinatesFromAbsoluteInoculum([n1,n2],args.AbsoluteCoordinates)))
                if (n1==0) and (n2==0):
                    break
            fp.write("\n")
        fp.close()

                
if __name__ == "__main__":
    main()
