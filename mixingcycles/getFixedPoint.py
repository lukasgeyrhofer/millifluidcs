#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.stats import poisson

def getGrowthMatrix(parameters,size = 300):
    def get_time_to_substrate_depletion(m1,m2,alpha=1.,maxsteps = 10000, precision = 1e-10):
        initialcells = np.array([m1,m2])
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
        
    mn1 = np.zeros((size,size))
    mn2 = np.zeros((size,size))
    for i in range(size):
        for j in range(size):
            t = get_time_to_substrate_depletion(i,j)
            mn1[i,j],mn2[i,j] = np.array([i,j]) * np.exp(p['growth']*t) * p['dilution']

    return mn1,mn2

def prob(m,n,cutoff):
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
parser_populations.add_argument("-Y","--yieldfactor",type=float,nargs="*",default=[1.,2.])
parser_populations.add_argument("-S","--substrateconcentration",type=float,default=1e4)
parser_populations.add_argument("-d","--dilutionfactor",type=float,default=2e-4)
parser_populations.add_argument("-T","--mixingtime",type=float,default=12.)

parser_algorithm = parser.add_argument_group(description = "=== Algorithm parameters ===")
parser_algorithm.add_argument("-N","--newtonraphson",action="store_true",default=False,help = "Plain iteration of dynamics or try to use NR to estimate fixed point")
parser_algorithm.add_argument("-m","--maxM",type=int,default=300,help = "maximum possible number for seeding[default: 300]")
parser_algorithm.add_argument("-M","--maxiterations",type=int,default=1000, help = "maximum number of iterations [default: 1000]")
parser_algorithm.add_argument("-A","--alpha",type=float,default=1.,help = "convergence parameter for NR [default: 1.0]")
parser_algorithm.add_argument("-p","--precision",type=float,default=1e-14,help = "relative precision as premature stopping condition, computed as sum( (dn/n)^2 ) [default: 1e-14]")
parser_algorithm.add_argument("-v","--verbose",action="store_true",default=False,help = "output current values every iteration step")
parser_algorithm.add_argument("-c","--cutoff",type=float,default=1e-100,help = "cutoff probabilities lower than this value [default: 1e-100]")
args = parser.parse_args()



# generate parameter dictionary
p = {}
p['yield']      = np.array(args.yieldfactor)
p['growth']     = np.array(args.growthrates)
p['substrate']  = args.substrateconcentration
p['mixingtime'] = args.mixingtime
p['dilution']   = args.dilutionfactor

# initialize necessary variables
growth1,growth2 = getGrowthMatrix(p,size = args.maxM)
m = np.arange(args.maxM)
j = np.zeros((2,2))

# initial condition are the respective fixed points on the axis
n = p['dilution']/(1-p['dilution'])*p['substrate']*p['yield']

for i in range(args.maxiterations):
    if args.verbose:
        print "#{:4d} {:12.8f} {:12.8f}".format(i,n[0],n[1])
    
    # probabilities for seeding new droplets, assumed to be poissonian
    px,py = prob(m,n,args.cutoff)
    
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
    
    # stop condition met? is the change small enough?
    if np.sum((dn[n>0]/n[n>0])**2) < args.precision:
        break

# final output
print "{:.6f} {:.6f} {:14.8f} {:14.8f} {:4d}".format(p['growth'][0]/p['growth'][1],p['yield'][0]/p['yield'][1],n[0],n[1],i)


