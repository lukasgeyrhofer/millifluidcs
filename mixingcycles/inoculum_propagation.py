#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import cairo
import scipy.misc as spm
import PIL

import growthclasses as gc


def inoc_from_coords(x,newcoordinates = False):
    if newcoordinates:
        return np.array([x[0]*x[1],x[0]*(1.-x[1])],dtype=np.float64)
    else:
        return np.array([x[0],x[1]],dtype=np.float64)


def value_to_color(value):
    


def write_image(filename,data,parameters):
    #CairoImage = cairo.ImageSurface(cairo.FORMAT_ARGB32,data.shape[1] * parameters['pixelsperpoint'],data.shape[0] * parameters['pixelsperpoint'])
    #context    = cairo.Context(CairoImage)

    if parameters.get('logscale',False):
        data  = np.log(data)
        data -= min(data)
        data /= max(data)
    elif parameters.get('rescale',False):
        data -= min(data)
        data /= max(data)
        
    img = msc.toimage(data,high = 1., low = 0.)
    img.save(filename,"PNG")
    
    
def verbose(msg,v=True):
    if v:
        if msg[-1] != '\n':
            msg += '\n'
        sys.stdout.write(msg)


def SGM_seeding(inocdens):
    newdens = np.zeros(np.shape(inocdens))
    for i,n1 in enumerate(m1):
        for j,n2 in enumerate(m2):
            if inocdens[i,j] > 0:
                p1 = gc.PoissonSeedingVectors(m1,[n1])[0]
                p2 = gc.PoissonSeedingVectors(m2,[n2])[0]
                
                newdens += inocdens[i,j] * np.outer(p1,p2)
    return newdens


def SGM_growth(inocdens,dilution):
    newdens = np.zeros(np.shape(inocdens))
    for i,n1 in enumerate(m1):
        for j,n2 in enumerate(m2):
            if inocdens[i,j] > 0:
                g1 = gm1[i,j] * dilution
                g2 = gm2[i,j] * dilution
                
                if g1 < m1[-1] and g2 < m2[-1]:
                    i1 = int(np.floor(g1))
                    i2 = int(np.floor(g2))
                    r1 = i1-g1
                    r2 = i2-g2
                    newdens[i1  ,i2  ] += inocdens[i,j] * (1-r1) * (1-r2)
                    newdens[i1+1,i2  ] += inocdens[i,j] *    r1  * (1-r2)
                    newdens[i1  ,i2+1] += inocdens[i,j] * (1-r1) *    r2
                    newdens[i1+1,i2+1] += inocdens[i,j] *    r1  *    r2
    return newdens


def SGM_mixing(inocdens,mode):
    if mode == "mixing":
        newdens = np.zeros(np.shape(inocdens))
        
        avg1 = np.dot(m1,np.sum(inocdens,axis=1))
        avg2 = np.dot(m2,np.sum(inocdens,axis=0))
        
        if avg1 < m1[-1] and avg2 < m2[-1]:
            i1 = int(np.floor(avg1))
            i2 = int(np.floor(avg2))
            r1 = i1-avg1
            r2 = i2-avg2
            
            newdens[i1  ,i2  ] += (1-r1) * (1-r2)
            newdens[i1+1,i2  ] +=    r1  * (1-r2)
            newdens[i1  ,i2+1] += (1-r1) *    r2
            newdens[i1+1,i2+1] +=    r1  *    r2

        return newdens
    else:
        return inocdens
                    

def __main__():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(descripton = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infile")
    parser_io.add_argument("-o","--outbasename",default="out")
    parser_io.add_argument("-v","--verbose",default=False,action = "store_true")

    parser_img = parser.add_argument_group(description = "==== Output image parameters ====")
    parser_img.add_argument("-p","--pixelsperpoint",type=int,default=2)

    parser_dilution = parser.add_argument_group(description = "==== Dilution values ====")
    parser_dilution.add_argument("-d","--dilutionmin",type=float,default=1e-6)
    parser_dilution.add_argument("-D","--dilutionmax",type=float,default=None)
    parser_dilution.add_argument("-K","--dilutionsteps",type=int,default=10)
    parser_dilution.add_argument("-L","--dilutionlogscale",default = False, action = "store_true")

    parser_lattice = parser.add_argument_group(description = "==== Lattice for cycles ====")
    parser_lattice.add_argument("-C","--newcoordinates",default=False,action="store_true",help="Use (n,x) instead of (n1,n2) as coordinates")
    parser_lattice.add_argument("-N","--maxInoculum",type=float,default=40)
    parser_lattice.add_argument("-n","--stepInoculum",type=float,default=2)
    parser_lattice.add_argument("-x","--stepFraction",type=float,default=.05)

    parser_mode = parser.add_argument_group(description = "==== Mode ====")
    parser_mode.add_argument("-m","--mode",choices = ("mixing","singlepop"),default = "mixing")
    parser_mode.add_argument("-M","--steps",default=100,type=int)
    parser_mode.add_argument("-I","--inoculum",default=[10,10],type=float,nargs=2)

    args=parser.parse_args()

    try:
        g = pickle.load(open(args.infile))
    except:
        raise IOError,"Could not open and load from pickle file"

    if not g.hasGrowthMatrix():
        raise ValueError,"Loaded pickle instance does not contain growthmatrix"

    verbose(g.ParameterString())

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
        xlist = nlist
        shape = (len(nlist),len(nlist))

    global gm1
    global gm2
    gm1   = g.growthmatrix[:,:,0]
    gm2   = g.growthmatrix[:,:,1]
    global m1
    global m2
    m1,m2 = g.growthmatrixgrid


    for dilution in dlist:
        # initialize with exact inoculum with single value (nearest the values provided as cmdline parameter)
        inocdens = np.zeros((len(m1),len(m2)))
        inocdens[np.argmin(m1[m1 >= inoc_from_coords(args.inoculum)[0]]),np.argmin(m2[m2 >= inoc_from_coords(args.inoculum)[1]])] = 1
        
        # write this as zeroth image
        outfn = args.outbasename + "_D{:.3e}_{%04d}".format(dilution,0)
        write_image(outfn,inocdens)
        
        for step in range(1,args.steps+1):
            
            verbose("D = {:.3e}, step {:4}".format(dilution,step))
            
            # update seeding probabilities from last cycle
            inocdens = SGM_seeding(inocdens)
            
            # store this seeding as image
            outfn = args.outbasename + "_D{:.3e}_{%04d}".format(dilution,step)
            write_image(outfn,inocdens)

            # growth
            inocdens = SGM_growth(inocdens,dilution)
            
            # dilution
            inocdens = SGM_mixing(inocdens,args.mode)
        
        


if __name__ == "__main__":
    __main__()
