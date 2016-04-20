#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.misc import factorial
from scipy.stats import poisson

# Return number of cells in droplet after growth and dilution for next round
def get_growth(initialcells,parameters):
    # Estimate time when resources are used up. This is an important quantity when mixtures of strains coexist in droplet.
    # Usually the fast growth strain uses up a large share of the whole amount of resources
    def get_time_to_substrate_depletion(alpha=1.,maxsteps = 10000, precision = 1e-10):
        # internal function to have all parameters.
        t0 = 0
        if np.sum(initialcells) > 0:
            # initial condition for iteration is assuming only strain with highest expected yield is present
            p = (1.*initialcells/parameters['yield']).argmax() # get ID of strain
            t1 = np.log(parameters['substrate']*parameters['yield'][p]/initialcells[p]+1.)/parameters['growth'][p] # set first estimate of time for single strain
            i = 0
            while ((t1-t0)/t1)**2 > precision**2:
                t0 = t1
                # Newton-Raphson iteration to refine solution
                t1 += alpha*(parameters['substrate']-np.sum(initialcells/parameters['yield']*(np.exp(parameters['growth']*t1)-1.)))/(np.sum(initialcells/parameters['yield']*parameters['growth']*np.exp(parameters['growth']*t1)))
                i+=1
                # should not iterate infinitely
                if i > maxsteps:
                    raise ValueError
            return min(t1,parameters['mixingtime'])
        else:
            return 0

    t = get_time_to_substrate_depletion()
    return initialcells * np.exp(parameters['growth'] * t) * parameters['dilution'],t



parser = argparse.ArgumentParser()
parser.add_argument("-m","--maxsize",type=int,default=300)
parser.add_argument("-n","--outputmax",type=float,default=None)
parser.add_argument("-D","--outputdx",type=float,default=1.)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-O","--poissonoutfile",default=None)
parser.add_argument("-T","--mixingtime",default=10,type=float)
parser.add_argument("-S","--substrateconcentration",default=1e4,type=float)
parser.add_argument("-d","--dilutionfactor",default=2e-4,type=float)
parser.add_argument("-a","--growthrates",type=float,nargs=2,default=[2.,1.])
parser.add_argument("-Y","--yieldfactor",type=float,nargs=2,default=[1.,2.])
args = parser.parse_args()

# generate parameter dictionary
p = {}
p['yield']      = np.array(args.yieldfactor)
p['growth']     = np.array(args.growthrates)
p['substrate']  = args.substrateconcentration
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilutionfactor


print p

# assume index [i,j] are numbers of cells in droplet:
# what is the expected final number of cells for all combinations of initial conditions?
# (up to args.maxsize cells of each type)
n0 = np.zeros((args.maxsize,args.maxsize,2))
t0 = np.zeros((args.maxsize,args.maxsize))
for i in range(args.maxsize):
    for j in range(args.maxsize):
        initialcells = 1.*np.array([i,j])
        n0[i,j,:],t0[i,j] = get_growth(initialcells,p)

# write outputfile
if args.outfile != None:
    fp = open(args.outfile,"w")
    for i in range(args.maxsize):
        for j in range(args.maxsize):
            print >> fp,"{:3d} {:3d} {:.5e} {:.5e} {:.5e}".format(i,j,n0[i,j,0],n0[i,j,1],t0[i,j])
        print >> fp
    fp.close()
elif args.poissonoutfile == None:
    # output to stdout if no outfile is specified at all
    for i in range(args.maxsize):
        for j in range(args.maxsize):
            print "{:3d} {:3d} {:.5e} {:.5e} {:.5e}".format(i,j,n0[i,j,0],n0[i,j,1],t0[i,j])
        print

# add noise:
# number of initial cells in droplet is drawn from poisson distribution
# thus, average distribution of growth with poissonian weights
if args.poissonoutfile != None:
    # get coordinates of initial conditions
    x = np.arange(0,args.outputmax if args.outputmax!=None else args.maxsize,args.outputdx)

    # first compute all relevant poisson distributions and store them in an array
    px = np.zeros((len(x),args.maxsize))
    px[0,0] = 1.
    for i in range(len(x)):
        px[i] = poisson.pmf(np.arange(args.maxsize),x[i])
    
    n1 = np.zeros((len(x),len(x),2))
    t1 = np.zeros((len(x),len(x)))

    fp = open(args.poissonoutfile,"w")
    for i in range(len(x)):
        for j in range(len(x)):
            # average values from above with poissonian weights
            n1[i,j] = np.dot(px[i],np.dot(px[j],n0))
            t1[i,j] = np.dot(px[i],np.dot(px[j],t0))
            print >> fp,"{:.5e} {:.5e} {:.5e} {:.5e} {:.5e}".format(x[i],x[j],n1[i,j,0],n1[i,j,1],t1[i,j])
        print >> fp
    fp.close()













