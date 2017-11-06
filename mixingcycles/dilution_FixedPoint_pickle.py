#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pickle

import growthclasses as gc

def re(x):
    return float(np.real(x))
def im(x):
    return float(np.imag(x))


parser = argparse.ArgumentParser()
#parser = gc.AddGrowthParameters(parser,dilution=True)
parser.add_argument("-i","--infile",help = "Pickle-file with stored parameters and growthmatrix")
parser.add_argument("-d","--dilution",type=float,default=1e-4,help = "dilution rate")

parser_algorithm = parser.add_argument_group(description = "=== Algorithm parameters ===")
parser_algorithm.add_argument("-N","--newtonraphson",action="store_true",default=False,help = "Plain iteration of dynamics or try to use NR to estimate fixed point")
#parser_algorithm.add_argument("-m","--maxM",type=int,default=100,help = "maximum possible number for seeding [default: 100]")
parser_algorithm.add_argument("-p","--precision",type=float,default=1e-20,help = "relative precision as premature stopping condition, computed as sum( (dn/n)^2 ) [default: 1e-20]")
parser_algorithm.add_argument("-M","--maxiterations",type=int,default=None, help = "maximum number of iterations [default: None, iterate until precision is reached]")
parser_algorithm.add_argument("-A","--alpha",type=float,default=1.,help = "convergence parameter for NR [default: 1.0]")
parser_algorithm.add_argument("-c","--cutoff",type=float,default=1e-100,help = "cutoff probabilities lower than this value [default: 1e-100]")

parser_general = parser.add_argument_group(description = "=== General and I/O parameters ===")
parser_general.add_argument("-v","--verbose",action="store_true",default=False,help = "output current values every iteration step")
parser_general.add_argument("-V","--printeigenvectors",default=False,action="store_true",help = "print eigenvectors of linearized iteration map")
parser_general.add_argument("-I","--initialconditions",default=None,nargs="*",help="Override initial conditions when set")
parser_general.add_argument("-C","--complexOutput",default=False,action="store_true",help="Print real and imaginary parts of eigenvalues (and eigenvectors) [default: only real]")
args = parser.parse_args()


try:
    g = pickle.load(open(args.infile))
except:
    raise IOError("could not open pickle file")


gm1   = g.growthmatrix[:,:,0]*args.dilution
gm2   = g.growthmatrix[:,:,1]*args.dilution
mx,my = g.growthmatrixgrid

# initial condition are the respective fixed points on the axis
if args.initialconditions is None:
    g.setDilution(args.dilution)
    n = g.getSingleStrainFixedPointsApproximate()
    
else:
    n   = np.array(args.initialconditions,dtype=float)
    assert len(n) == g.numstrains

dn         = n                  # increase in one step
j          = np.zeros((2,2))    # jacobian
stepcount  = 0                  # number of steps for debugging

if args.verbose:
    sys.stdout.write("# starting iterations ... \n")

while np.sum((dn[n>0]/n[n>0])**2) > args.precision:
    if args.verbose:
        sys.stderr.write("{:4d} {:12.8e} {:12.8e}\n".format(stepcount,n[0],n[1]))
    
    # probabilities for seeding new droplets, assumed to be poissonian
    px,dpx = gc.PoissonSeedingVectors(mx,[n[0]],cutoff = args.cutoff,diff=True)
    py,dpy = gc.PoissonSeedingVectors(my,[n[1]],cutoff = args.cutoff,diff=True)
    
    # construct iteration function for growth and dilution
    # by weighting growth with the probability of how droplets are seeded
    growth1 = np.dot(py[0], np.dot(px[0], gm1))
    growth2 = np.dot(py[0], np.dot(px[0], gm2))
    fn = np.array([growth1,growth2]) - n
    
    if args.newtonraphson:
        # NR iterations 

        # get jacobian of dynamics
        j[0,0] = np.dot( py[0], np.dot(dpx[0], gm1)) - 1.
        j[0,1] = np.dot(dpy[0], np.dot( px[0], gm1))
        j[1,0] = np.dot( py[0], np.dot(dpx[0], gm2))
        j[1,1] = np.dot(dpy[0], np.dot( px[0], gm2)) - 1.
        
        # calculate step in NR iteration
        dn = -args.alpha * np.dot(np.linalg.inv(j),fn)
    
    else:
        # simple iteration of the function, hoping it converges at some point
        dn = fn
    
    # apply changes
    n += dn
    n[n<0] = 0
    
    if not args.maxiterations is None:
        if i > args.maxiterations:
            break
    stepcount += 1


# stability of fixed point is checked with jacobian
px,dpx = gc.PoissonSeedingVectors(mx,[n[0]],cutoff = args.cutoff,diff=True)
py,dpy = gc.PoissonSeedingVectors(my,[n[1]],cutoff = args.cutoff,diff=True)

j[0,0] = np.dot( py[0], np.dot(dpx[0], gm1))
j[0,1] = np.dot(dpy[0], np.dot( px[0], gm1))
j[1,0] = np.dot( py[0], np.dot(dpx[0], gm2))
j[1,1] = np.dot(dpy[0], np.dot( px[0], gm2))


w,v = np.linalg.eig(j)

# final output

outputstring  = "{:14.6f} {:14.6f} {:14.6e} {:14.6e} {:4d}".format(g.growthrates[1]/g.growthrates[0], g.yieldfactors[1]/g.yieldfactors[0], n[0], n[1], stepcount)
outputstring += " {:14.6f} {:14.6f}".format(re(w[0]),re(w[1]))
if args.complexOutput:
    # have yet to find complex eigenvalues
    outputstring += " {:14.6f} {:14.6f}".format(im(w[0]),im(w[1]))

if args.printeigenvectors:
    outputstring += " {:14.6f} {:14.6f} {:14.6f} {:14.6f}".format(re(v[0,0]),re(v[1,0]),re(v[0,1]),re(v[1,1]))
    if args.complexOutput:
        outputstring += " {:11.6f} {:11.6f} {:11.6f} {:11.6f}".format(im(v[0][0]),im(v[0][1]),im(v[1][0]),im(v[1][1]))
sys.stdout.write(outputstring + "\n")


