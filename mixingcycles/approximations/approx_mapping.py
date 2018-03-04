#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import itertools
from skimage import measure

sys.path.append(sys.path[0] + '/..')
import growthclasses as gc



def verbose(msg,verbose = True):
    if verbose:
        if msg[-1] != '\n':
            msg += '\n'
        sys.stdout.write(msg)


def correction_term(m1,m2,model,modelparameters):
    x = 0
    if m1+m2>0:x=float(m1)/(m1+m2)
    r = 1 + modelparameters['dy'][0] * (2.*x-1.)
    if model == 'AB':
        if m1 * modelparameters['cmdline'][0] + m2 * modelparameters['cmdline'][1] < 1:
            r = 0
    elif model == 'PVD':
        r *= modelparameters['cmdline'][0]
    return r


def coord_to_inoc(c):
    return np.array([c[1]*c[0],(1.-c[1])*c[0]])


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
parser = gc.AddGrowthParameters(parser)

parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-o","--baseoutfilename",                        default = "out")
parser_io.add_argument("-v","--verbose",         action = "store_true", default = False)
parser_io.add_argument("-C","--computationmode", choices=['SS','NC'],   default = 'SS') # 'SS': single step dynamics ,'NC': nullclines

parser_dilution = parser.add_argument_group(description = "==== Parameters for dilution values ====")
parser_dilution.add_argument("-d","--dilutionmin",   type=float, default = 1.e-4)
parser_dilution.add_argument("-D","--dilutionmax",   type=float, default = None)
parser_dilution.add_argument("-K","--dilutionsteps", type=int,   default = 10)
parser_dilution.add_argument("-L","--dilutionlogscale",          default = False, action = "store_true")

parser_lattice = parser.add_argument_group(description = "==== Output lattice ====")
parser_lattice.add_argument("-m","--maxM",   type=int,              default = 100)
parser_lattice.add_argument("-n","--minN",   type=float,            default = 1.e-2)
parser_lattice.add_argument("-N","--maxN",   type=float,            default = 50.)
parser_lattice.add_argument("-k","--stepsN", type=int,              default = 101)
parser_lattice.add_argument("-l","--logN",   action = "store_true", default = False)
parser_lattice.add_argument("-x","--stepsX", type=int,              default = 21)

parser_model = parser.add_argument_group(description = "==== Within droplet dynamics ====")
parser_model.add_argument("-M","--model",           choices = ['GY','PVD','AB'], default = 'GY')
parser_model.add_argument("-p","--modelparameters", nargs = "*", type = float,   default = [])

args = parser.parse_args()



verbose("# generating initial conditions",args.verbose)
g    = gc.GrowthDynamics(**vars(args))

if args.logN:   nlist = np.exp(np.linspace(start = np.log(args.minN),stop = np.log(args.maxN), num = args.stepsN))
else:           nlist = np.linspace(start = args.minN, stop = args.maxN, num = args.stepsN)
xlist                 = np.linspace(start = 0, stop = 1, num  = args.stepsX)
mlist                 = np.arange  (start = 0, stop = args.maxM, dtype=int)

nmat                  = np.repeat(np.expand_dims(nlist, axis = 1), repeats = len(xlist), axis = 1)
xmat                  = np.repeat(np.expand_dims(xlist, axis = 0), repeats = len(nlist), axis = 0)

gmshape               = (args.maxM,args.maxM)
outshape              = (len(nlist),len(xlist))

f1                    = np.zeros(shape = gmshape)
f2                    = np.zeros(shape = gmshape)

if args.dilutionmax is None:
    dlist = np.array([args.dilutionmin])
else:
    if args.dilutionlogscale:   dlist = np.exp(np.linspace(start = np.log(args.dilutionmin), stop = np.log(args.dilutionmax), num = args.dilutionsteps))
    else:                       dlist = np.linspace(start = args.dilutionmin, stop = args.dilutionmax, num = args.dilutionsteps)

y     = np.mean(g.yieldfactors)
a     = np.mean(g.growthrates)
dy    = (g.yieldfactors - y)/y
da    = (g.growthrates  - a)/a
syda  = np.power(g.env.substrate * y,da[0])

modelparameters = {'cmdline':args.modelparameters,'dy':dy,'da':da}

compute_only_singlestep = (args.computationmode == 'SS')



verbose("# generating approximations for growth terms",args.verbose)
for m1,m2 in itertools.product(mlist,repeat=2):
    n = float(m1+m2)
    if n > 0:
        x = m1/n
        f1[m1,m2] = x     * np.power(n,-da[0]) * np.power(correction_term(m1,m2,args.model,modelparameters),1+da[0])
        f2[m1,m2] = (1-x) * np.power(n,+da[0]) * np.power(correction_term(m1,m2,args.model,modelparameters),1-da[0])



verbose("# iterating dilutions",args.verbose)
for c,dilution in enumerate(dlist):
    verbose("# D = {:.3e}\n".format(dilution),args.verbose)

    if compute_only_singlestep:
        fp = open(args.baseoutfilename + '_D{:.3e}'.format(dilution),"w")

    newn   = np.zeros(shape = outshape)
    newx   = np.zeros(shape = outshape)
    avg_f1 = np.zeros(shape = outshape)
    avg_f2 = np.zeros(shape = outshape)
    
    for i,n in enumerate(nlist):
        for j,x in enumerate(xlist):
            p = gc.PoissonSeedingVectors(mlist,coord_to_inoc(np.array([n,x])))

            avg_f1[i,j] = np.dot(p[1],np.dot(p[0],f1))
            avg_f2[i,j] = np.dot(p[1],np.dot(p[0],f2))
            
            if avg_f1[i,j] + avg_f2[i,j]/(syda*syda) > 0:
                newn[i,j] = dilution * g.env.substrate * y * syda * (avg_f1[i,j] + avg_f2[i,j]/(syda*syda))
                if c == 0:
                    # only need to compute for first dilution value, this does not change thereafter
                    newx[i,j] = avg_f1[i,j]/(avg_f1[i,j] + avg_f2[i,j]/(syda*syda))
            
            if compute_only_singlestep: fp.write("{:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n".format(coord[0],coord[1],newn[i,j],newx[i,j],avg_f1[i,j],avg_f2[i,j]))
        if compute_only_singlestep:     fp.write("\n")
    if compute_only_singlestep:         fp.close()



    if not compute_only_singlestep:
        # get also nullclines
        verbose("# single step computations finished, computing nullclines",args.verbose)
        
        if c == 0:
            # only need to compute for first dilution value, this does not change thereafter
            cont_x = measure.find_contours(newx - xmat,0)
            write_contours_to_file(cont_x,args.baseoutfilename + '_X',nlist,xlist)
        
        cont_n = measure.find_contours(newn - nmat,0)
        write_contours_to_file(cont_n,args.baseoutfilename + '_N_D{:.3e}'.format(dilution),nlist,xlist)
        
    


