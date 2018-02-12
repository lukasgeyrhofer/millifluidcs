#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle
from skimage import measure


import growthclasses as gc


def write_contours_to_file(contours, filename, nlist, xlist = None):
    fp = open(filename,"w")
    for c in contours:
        for i in range(len(c)):
            ix = int(np.floor(c[i,0]))
            iy = int(np.floor(c[i,1]))
            px = c[i,0] - ix
            py = c[i,1] - iy
            try:        cx = (1.-px)*nlist[ix] + px*nlist[ix+1]
            except:     cx = nlist[ix]
            if xlist is None:
                try:    cy = (1.-py)*nlist[iy] + py*nlist[iy+1]
                except: cy = nlist[iy]
            else:
                try:    cy = (1.-py)*xlist[iy] + py*xlist[iy+1]
                except: cy = xlist[iy]
            fp.write('{:14.6e} {:14.6e}\n'.format(cx,cy))
        fp.write('\n')
    fp.close()


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


parser_flowmap = parser.add_argument_group(description = "==== Flowmap between mixing cycles ====")
parser_flowmap.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
parser_flowmap.add_argument("-N","--maxInoculum",type=float,default=40)
parser_flowmap.add_argument("-n","--stepInoculum",type=float,default=2)
parser_flowmap.add_argument("-x","--stepFraction",type=float,default=.05)

args=parser.parse_args()

try:
    g = pickle.load(open(args.infile))
except:
    raise IOError,"Could not open and load from pickle file"

if not g.hasGrowthMatrix():
    raise ValueError,"Loaded pickle instance does not contain growthmatrix"

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

if args.newcoordinates:
    nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    xlist = np.arange(start = 0,stop = 1 + .5*args.stepFraction,step = args.stepFraction)
    shape = (len(nlist),len(xlist))
else:
    nlist = np.arange(start = 0,stop = args.maxInoculum,step = args.stepInoculum)
    xlist = None
    shape = (len(nlist),len(nlist))


m1,m2 = g.growthmatrixgrid
gm1   = g.growthmatrix[:,:,0]
gm2   = g.growthmatrix[:,:,1]


g1  = np.zeros(shape,dtype=np.float64)
g2  = np.zeros(shape,dtype=np.float64)
rr1 = np.zeros(shape,dtype=np.float64)
r1  = np.zeros(shape,dtype=np.float64)

if not args.newcoordinates:
    for i,n1 in enumerate(nlist):
        for j,n2 in enumerate(nlist):
            p1 = gc.PoissonSeedingVectors(m1,[n1])
            p2 = gc.PoissonSeedingVectors(m2,[n2])
            g1[i,j] = np.dot(p2[0],np.dot(p1[0],gm1))
            g2[i,j] = np.dot(p2[0],np.dot(p1[0],gm2))

    sn1 = np.array(np.repeat(np.transpose([nlist]),len(nlist),axis=1),dtype=float)
    sn2 = np.repeat(np.array([nlist],dtype = float),len(nlist),axis=0)

    rr1[g1+g2>0]  = (g1[g1+g2>0])/((g1+g2)[g1+g2>0])
    r1[sn1+sn2>0] = (sn1[sn1+sn2>0])/((sn1+sn2)[sn1+sn2>0])

else:
    for i,n in enumerate(nlist):
        for j,x in enumerate(xlist):
            p1 = gc.PoissonSeedingVectors(m1,[x*n])
            p2 = gc.PoissonSeedingVectors(m2,[(1-x)*n])
            g1[i,j] = np.dot(p2[0],np.dot(p1[0],gm1))
            g2[i,j] = np.dot(p2[0],np.dot(p1[0],gm2))
            
    sn1 = np.outer(nlist,xlist)
    sn2 = np.outer(nlist,1-xlist)

    rr1[g1+g2>0] = (g1[g1+g2>0])/((g1+g2)[g1+g2>0])
    r1[sn1+sn2>0] = (sn1[sn1+sn2>0])/((sn1+sn2)[sn1+sn2>0])


if args.verbose:
    sys.stdout.write('\n computing nullcline for fraction of strains\n')
cont_xx = measure.find_contours(rr1 - r1,0)
write_contours_to_file(cont_xx,args.baseoutfilename + '_X',nlist,xlist)

for dilution in dlist:
    if args.verbose:
        sys.stdout.write(' computing nullclines for dilution D = {:.4e}\n'.format(dilution))
    cont_nn = measure.find_contours((g1 + g2) * dilution - sn1 - sn2,0)
    write_contours_to_file(cont_nn,args.baseoutfilename + '_N_D{:.3e}'.format(dilution),nlist,xlist)
    if args.OutputSinglestrainNullclines:
        cont_n1 = measure.find_contours(g1 * dilution - sn1,0)
        cont_n2 = measure.find_contours(g2 * dilution - sn2,0)
        write_contours_to_file(cont_n1,args.baseoutfilename + '_1_D{:.3e}'.format(dilution),nlist,xlist)
        write_contours_to_file(cont_n2,args.baseoutfilename + '_2_D{:.3e}'.format(dilution),nlist,xlist)

