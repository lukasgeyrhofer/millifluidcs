#!/usr/bin/env python3
#-*- coding: utf-8 -*-

'''
==================================
=  mixingcycles_Isoclines.py
==================================

    Computes isoclines of the dynamics of average inoculum sizes over multiple cycles


    Lukas Geyrhofer, l.geyrhofer@technion.ac.il, 2018

'''


import numpy as np
import argparse
import sys,math
import pickle
from skimage import measure
import growthclasses as gc


def write_contours_to_file(contours, filename, axis1, axis2):
    fp = open(filename,"w")
    for c in contours:
        for i in range(len(c)):
            ix = int(np.floor(c[i,0]))
            iy = int(np.floor(c[i,1]))
            px = c[i,0] - ix
            py = c[i,1] - iy
            try:    cx = (1.-px)*axis1[ix] + px*axis1[ix+1]
            except: cx = axis1[ix]
            try:    cy = (1.-py)*axis2[iy] + py*axis2[iy+1]
            except: cy = axis2[iy]
            fp.write('{:14.6e} {:14.6e}\n'.format(cx,cy))
        fp.write('\n')
    fp.close()


def main():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O ====")
    parser_io.add_argument("-i","--infile",required=True)
    parser_io.add_argument("-o","--baseoutfilename",default="out")
    parser_io.add_argument("-v","--verbose",action="store_true",default=False)
    parser_io.add_argument("-S","--OutputSinglestrainNullclines",action="store_true",default=False)

    parser_dilution = parser.add_argument_group(description = "==== Dilution values ====")
    parser_dilution.add_argument("-d","--dilutionmin",type=float,default=1e-6)
    parser_dilution.add_argument("-D","--dilutionmax",type=float,default=None)
    parser_dilution.add_argument("-K","--dilutionsteps",type=int,default=10)
    parser_dilution.add_argument("-L","--dilutionlogscale",default = False, action = "store_true")


    parser_lattice = parser.add_argument_group(description = "==== Lattice parameters ====")
    parser_lattice.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
    parser_lattice.add_argument("-N","--maxInoculum",type=float,default=40)
    parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
    parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)

    args=parser.parse_args()

    try:
        g = pickle.load(open(args.infile,'rb'), encoding = 'bytes')
    except:
        raise IOError("Could not open and load from pickle file")

    if not g.hasGrowthMatrix():
        raise ValueError("Loaded pickle instance does not contain growthmatrix")

    if args.verbose:
        sys.stdout.write(g.ParameterString())
        sys.stdout.write("\n generating matrices\n")

    if args.dilutionmax is None:
        dlist = np.array([args.dilutionmin])
    else:
        if args.dilutionlogscale:
            dlist = np.power(10,np.linspace(start = np.log10(args.dilutionmin),stop = np.log10(args.dilutionmax), num = args.dilutionsteps))
        else:
            dlist = np.linspace(start = args.dilutionmin,stop = args.dilutionmax,num = args.dilutionsteps)

    # get new axes, which depends on parameters above (in lattice parameter group)
    axis1,axis2 = gc.getInoculumAxes(**vars(args)) # either (n1,n2) or [ (n,x) if args.newcoordinates == True ]
    shape       = (len(axis1),len(axis2))

    # loaded from pickle file
    m1,m2       = g.growthmatrixgrid
    gm1         = g.growthmatrix[:,:,0]
    gm2         = g.growthmatrix[:,:,1]

    # matrices to store averages
    g1          = np.zeros(shape,dtype=np.float64) # avg'd growth strain 1
    g2          = np.zeros(shape,dtype=np.float64) # avg'd growth strain 2
    rr1         = np.zeros(shape,dtype=np.float64) # avg'd ratio of strains at end
    r1          = np.zeros(shape,dtype=np.float64) # avg'd ratio of strains at beginning
    sn1         = np.zeros(shape,dtype=np.float64) # number of cells of strain 1 in new matrix shape
    sn2         = np.zeros(shape,dtype=np.float64) # number of cells of strain 1 in new matrix shape

    # get all averages and store them in the appropriate matrices
    for i,a1 in enumerate(axis1):
        for j,a2 in enumerate(axis2):
            sn1[i,j],sn2[i,j] = gc.getAbsoluteInoculumNumbers([a1,a2],args.newcoordinates)
            g1[i,j] = gc.SeedingAverage(gm1, [sn1[i,j],sn2[i,j]])
            g2[i,j] = gc.SeedingAverage(gm2, [sn1[i,j],sn2[i,j]])

    rr1[g1+g2>0]  = (g1[g1+g2>0])/((g1+g2)[g1+g2>0])
    r1[sn1+sn2>0] = (sn1[sn1+sn2>0])/((sn1+sn2)[sn1+sn2>0])

    # output
    if args.verbose:
        sys.stdout.write('\n computing nullcline for fraction of strains\n')
    cont_xx = measure.find_contours(rr1 - r1,0)
    write_contours_to_file(cont_xx,args.baseoutfilename + '_X',axis1,axis2)

    for dilution in dlist:
        if args.verbose:
            sys.stdout.write(' computing nullclines for dilution D = {:.4e}\n'.format(dilution))
        cont_nn = measure.find_contours((g1 + g2) * dilution - sn1 - sn2,0)
        write_contours_to_file(cont_nn,args.baseoutfilename + '_N_D{:.3e}'.format(dilution),axis1,axis2)
        if args.OutputSinglestrainNullclines:
            cont_n1 = measure.find_contours(g1 * dilution - sn1,0)
            cont_n2 = measure.find_contours(g2 * dilution - sn2,0)
            write_contours_to_file(cont_n1,args.baseoutfilename + '_1_D{:.3e}'.format(dilution),axis1,axis2)
            write_contours_to_file(cont_n2,args.baseoutfilename + '_2_D{:.3e}'.format(dilution),axis1,axis2)

                
if __name__ == "__main__":
    main()
