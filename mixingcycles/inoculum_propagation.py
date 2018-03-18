#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import cairo
import scipy.misc as spm
import PIL
import pickle

import growthclasses as gc


def inoc_from_coords(x,newcoordinates = False):
    if newcoordinates:
        return np.array([x[0]*x[1],x[0]*(1.-x[1])],dtype=np.float64)
    else:
        return np.array([x[0],x[1]],dtype=np.float64)


def write_image(filename,data,parameters = dict()):
    #CairoImage = cairo.ImageSurface(cairo.FORMAT_ARGB32,data.shape[1] * parameters['pixelsperpoint'],data.shape[0] * parameters['pixelsperpoint'])
    #context    = cairo.Context(CairoImage)
    
    tmp = np.copy(data)

    if parameters.get('logscale',False):
        tmp  = np.log(tmp)
        if not parameters.get('logmin',None) is None:
            tmp[tmp < parameters.get('logmin',None)] = parameters.get('logmin',None)
        tmp[np.isnan(tmp)] = np.nanmin(tmp)
        tmp /= np.nanmin(tmp)
        tmp = 1 - tmp
    
    if parameters.get('rescale',False):
        tmp[np.isnan(tmp)] = np.nanmin(tmp)
        tmp -= np.nanmin(tmp)
        m    = np.nanmax(tmp)
        if m > 0: tmp /= m
    
    img = PIL.Image.fromarray((255 * tmp).astype('uint8'))
    img.save(filename)
    
    
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
    # check the growth of each inoculum, and add it to the seeding probabilities for the next cycle
    newdens = np.zeros(np.shape(inocdens))

    for i,n1 in enumerate(m1):
        for j,n2 in enumerate(m2):
            if inocdens[i,j] > 0:
                growth = np.array([gm1[i,j],gm2[i,j]],dtype=np.float64) * dilution
                newdens += inocdens[i,j] * inocdens_setsinglevalue(growth)
                
    return newdens


def SGM_mixing(inocdens,mode):
    if mode == "mixing":
        # compute average over all possible inocula, and set the whole inoculum density to this average
        avg = avg_inocdens(inocdens)
        return inocdens_setsinglevalue(avg)
    elif mode == "singlepop":
        # do not change inoculum distribution and let it evolve freely
        return np.copy(inocdens)
    else:
        return None


def inocdens_setsinglevalue(singleinoc,newcoordinates = False):
    # distribute single value of an averaged inoculum to the 4 most adjacent bins
    newdens = np.zeros((len(m1),len(m2)))
    ic = inoc_from_coords(singleinoc,newcoordinates)

    if ic[0] < m1[-1] and ic[1] < m2[-1]:
        i1 = len(m1[m1 <= ic[0]])-1
        i2 = len(m2[m2 <= ic[1]])-1
        r1 = (ic[0] - m1[i1])/(m1[i1+1] - m1[i1])
        r2 = (ic[1] - m2[i2])/(m2[i2+1] - m2[i2])
        
        newdens[i1  ,i2  ] += (1-r1) * (1-r2)
        newdens[i1+1,i2  ] +=    r1  * (1-r2)
        newdens[i1  ,i2+1] += (1-r1) *    r2
        newdens[i1+1,i2+1] +=    r1  *    r2

    return newdens


def avg_inocdens(inocdens):
    avg1 = np.dot(m1,np.sum(inocdens,axis=1))
    avg2 = np.dot(m2,np.sum(inocdens,axis=0))
    return np.array([avg1,avg2],dtype=np.float64)


def __main__():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infile")
    parser_io.add_argument("-o","--outbasename",default="out")
    parser_io.add_argument("-v","--verbose",default=False,action = "store_true")

    parser_img = parser.add_argument_group(description = "==== Output image parameters ====")
    parser_img.add_argument("-p","--pixelsperpoint",type=int,default=2) # unimplemented
    parser_img.add_argument("-l","--DensityLogscale",default=False,action="store_true")
    parser_img.add_argument("-P","--DensityLogscaleMin",default=None,type=float)
    parser_img.add_argument("-r","--DensityRescale",default=False,action="store_true")

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

    args = parser.parse_args()
    
    imageparameters = {
            'logscale': args.DensityLogscale,
            'rescale':  args.DensityRescale,
            'logmin':   args.DensityLogscaleMin
            }
    

    try:
        g = pickle.load(open(args.infile))
    except:
        raise IOError,"Could not open and load from pickle file"

    if not g.hasGrowthMatrix():
        raise ValueError,"Loaded pickle instance does not contain growthmatrix"

    verbose(g.ParameterString() + '\n\n',args.verbose)

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
        inocdens = inocdens_setsinglevalue(args.inoculum,args.newcoordinates)
        
        avg = avg_inocdens(inocdens)
        verbose("D = {:.3e}, step {:4}, avg [{:10.3e}, {:10.3e}]".format(dilution,0,avg[0],avg[1]),args.verbose)
        

        # write this as zeroth image
        outfn = args.outbasename + "_D{:.3e}_{:04d}.png".format(dilution,0)
        write_image(outfn,inocdens,imageparameters)
        
        for step in range(1,args.steps+1):
            
            
            # update seeding probabilities from last cycle
            inocdens = SGM_seeding(inocdens)
            
            # store this seeding as image
            outfn = args.outbasename + "_D{:.3e}_{:04d}.png".format(dilution,step)
            write_image(outfn,inocdens,imageparameters)

            # growth
            inocdens = SGM_growth(inocdens,dilution)
            
            # dilution
            inocdens = SGM_mixing(inocdens,args.mode)
            
            avg = avg_inocdens(inocdens)
            verbose("D = {:.3e}, step {:4}, avg [{:10.3e}, {:10.3e}]".format(dilution,step,avg[0],avg[1]),args.verbose)
        
        


if __name__ == "__main__":
    __main__()
