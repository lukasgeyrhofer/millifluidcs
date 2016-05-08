#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

from growthclasses import growthdynamics


def prob(m,n,cutoff = 1e-100):
    if n[0] > 0:
        px = poisson.pmf(m,n[0])
        px[px<cutoff] = 0.
        px[-1] += (1. - np.sum(px))
    else:
        px = np.zeros(len(m))
        px[0] = 1
    if n[1] > 0:
        py = poisson.pmf(m,n[1])
        py[py<cutoff] = 0.
        py[-1] += (1. - np.sum(py))
    else:
        py = np.zeros(len(m))
        py[0] = 1
    return px,py


parser = argparse.ArgumentParser()
parser_populations = parser.add_argument_group(description = "=== Parameters for population growth and mixing ===")
parser_populations.add_argument("-a","--growthrates",type=float,nargs="*",default=[2.,1.])
parser_populations.add_argument("-Y","--yieldrates",type=float,nargs="*",default=[1.,2.])
parser_populations.add_argument("-S","--substrateconcentration",type=float,default=1e4)
parser_populations.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
parser_populations.add_argument("-T","--mixingtime",type=float,default=12.)

parser_algorithm = parser.add_argument_group(description = "=== Algorithm parameters ===")
parser_algorithm.add_argument("-N","--newtonraphson",action="store_true",default=False,help = "Plain iteration of dynamics or try to use NR to estimate fixed point")
parser_algorithm.add_argument("-m","--maxM",type=int,default=300,help = "maximum possible number for seeding [default: 300]")
parser_algorithm.add_argument("-p","--precision",type=float,default=1e-14,help = "relative precision as premature stopping condition, computed as sum( (dn/n)^2 ) [default: 1e-14]")
parser_algorithm.add_argument("-M","--maxiterations",type=int,default=None, help = "maximum number of iterations [default: None, iterate until precision is reached]")
parser_algorithm.add_argument("-A","--alpha",type=float,default=1.,help = "convergence parameter for NR [default: 1.0]")
parser_algorithm.add_argument("-v","--verbose",action="store_true",default=False,help = "output current values every iteration step")
parser_algorithm.add_argument("-c","--cutoff",type=float,default=1e-100,help = "cutoff probabilities lower than this value [default: 1e-100]")
args = parser.parse_args()


# initialize necessary variables
if args.verbose: print >> sys.stderr,"# initializing growth matrix ..."

g = growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilutionfactor, substrate = args.substrateconcentration)

#t = g.getTimeToDepletionMatrix(args.maxM)
#print t
#exit(1)

growth1,growth2 = g.getGrowthMatrix(size = args.maxM)
# initial condition are the respective fixed points on the axis
n = g.getSingleStrainFixedPoints()


m = np.arange(args.maxM)
j = np.zeros((2,2))
dn = n
i = 0

if args.verbose: print >> sys.stderr,"# starting iterations ..."
while np.sum((dn[n>0]/n[n>0])**2) > args.precision:
    if args.verbose:
        print >> sys.stderr,"{:4d} {:12.8e} {:12.8e}".format(i,n[0],n[1])
    
    # probabilities for seeding new droplets, assumed to be poissonian
    px,py = prob(m,n,cutoff = args.cutoff)
    
    # construct iteration function for growth and dilution
    # by weighting growth with the probability of how droplets are seeded
    fn = np.array([np.dot(py,np.dot(px,growth1)),np.dot(py,np.dot(px,growth2))])
    
    if args.newtonraphson:
        # NR iterations 
        # generate derivatives with respect to variables for NR, dP[m|n]/dn
        if n[0] > 0:    dpx = (m/n[0] - 1.)*px
        else:           dpx = -px
        if n[1] > 0:    dpy = (m/n[1] - 1.)*py
        else:           dpy = -py

        # get jacobian of dynamics
        j[0,0] = np.dot(dpy,np.dot(px,growth1))-1.
        j[0,1] = np.dot(py,np.dot(dpx,growth1))
        j[1,0] = np.dot(dpy,np.dot(px,growth2))
        j[1,1] = np.dot(py,np.dot(dpx,growth2))-1.
        
        # calculate step in NR iteration
        dn = -args.alpha * np.dot(np.linalg.inv(j),fn-n)
    
    else:
        # simple iteration of the function, hoping it converges at some point
        dn = fn-n
    
    # apply changes
    n += dn
    n[n<0] = 0
    
    if not args.maxiterations is None:
        if i<args.maxiterations:    i += 1
        else:                       break

# final output
print "{:.6f} {:.6f} {:14.8e} {:14.8e} {:4d}".format(g.growthrates[0]/g.growthrates[1],g.yieldrates[0]/g.yieldrates[1],n[0],n[1],i)


