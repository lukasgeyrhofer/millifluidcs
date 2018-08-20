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
    
    parser_flowmap = parser.add_argument_group(description = "==== Lattice parameters ====")
    parser_flowmap_startingconditions = parser_flowmap.add_mutually_exclusive_group()
    parser_flowmap_startingconditions.add_argument("-n","--maxIC",type=float,default=40)
    parser_flowmap_startingconditions.add_argument("-N","--initialcoordinatesfile",default=None)
    parser_flowmap.add_argument("-s","--stepIC",type=float,default=2)
    parser_flowmap.add_argument("-l","--trajectorylength",type=int,default=20)

    args = parser.parse_args()


    g     = gc.LoadGM(**vars(args))
    dlist = gc.getDilutionList(**vars(args))

    if args.initialcoordinatesfile is None:
        nlist = np.arange(start = 0,stop = args.maxIC,step = args.stepIC)
        coordinates = list(itertools.product(nlist,repeat=2))
    else:
        try:
            fp_coords = open(args.initialcoordinatesfile)
        except:
            raise IOError("could not open file '{}' to load coordinates".format(args.initialcoordinatesfile))
        coordinates = list()
        for line in fp_coords.readlines():
            try:
                values = np.array(line.split(),dtype=np.float)
                if len(values) >= 2:
                    coordinates.append(values[:2])
            except:
                continue
            fp_coords.close()


    mx,my = g.growthmatrixgrid
    gm1   = g.growthmatrix[:,:,0]
    gm2   = g.growthmatrix[:,:,1]

    if args.verbose:
        sys.stderr.write(str(g)

    for dilution in dlist:
        if args.verbose:
            sys.stderr.write("# computing trajectories for D = {:e}\n".format(dilution))
        fp = open(args.outfile + "_D{:.3e}".format(dilution),"w")
        for icx,icy in coordinates:
            x = icx
            y = icy
            fp.write("{} {}\n".format(x,y))
            for i in range(args.trajectorylength):
                px = gc.PoissonSeedingVectors(mx,x)[0]
                py = gc.PoissonSeedingVectors(my,y)[0]
                x = np.dot(py,np.dot(px,gm1))*dilution
                y = np.dot(py,np.dot(px,gm2))*dilution
                fp.write("{} {}\n".format(x,y))
                if (x==0) and (y==0):
                    break
            fp.write("\n")
        fp.close()

                
if __name__ == "__main__":
    main()
