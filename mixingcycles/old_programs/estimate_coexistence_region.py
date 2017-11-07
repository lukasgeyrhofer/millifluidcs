#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc

from scipy.special import lambertw

def lw(z):
    return float(np.real(lambertw(z)))

def gamma1(m1,m2,gr,yi,di,su):
    e1 = su/m1
    e2 = su*yi/m2
    return e2/e1/(gr-1)/lw(e2/(gr-1) * np.exp((1+e2/e1)/(gr-1)))


parser = argparse.ArgumentParser()

parser.add_argument("-a","--amin",type=float,default=.1)
parser.add_argument("-A","--amax",type=float,default=10)
parser.add_argument("-y","--ymin",type=float,default=.1)
parser.add_argument("-Y","--ymax",type=float,default=10)
parser.add_argument("-L","--logscale",default = False, action = "store_true")
parser.add_argument("-n","--gridsize",type=int,default=100)

parser.add_argument("-m","--maxM",type=int,default=30)

parser.add_argument("-d","--dilution",type=float, default=5e-6)
parser.add_argument("-S","--substrate",type=float,default=4e5)

args = parser.parse_args()

g = gc.GrowthDynamics(numstrains=2)
g.setDilution(args.dilution)
g.setSubstrate(args.substrate)
g.setMixingTime(100)

if args.logscale:
    alist = np.power(10,np.linspace(args.amin,args.amax,num = args.gridsize))
    ylist = np.power(10,np.linspace(args.ymin,args.ymax,num = args.gridsize))
else:
    alist = np.linspace(args.amin,args.amax,num = args.gridsize)
    ylist = np.linspace(args.ymin,args.ymax,num = args.gridsize)

dSY1 = args.dilution * args.substrate

for a in alist:
    g.growthrates = np.array([1,a])
    for y in ylist:
        g.yieldfactors = np.array([1,y])
        #gm = g.getGrowthMatrix(size = args.maxM)
        fp_exact = g.getSingleStrainFixedPointsPoissonSeeding()
        fp = args.dilution * args.substrate * np.array([1,y])
        g1_exact = g.Growth(initialcells = np.array([fp_exact[0],1]))[0]/g.Growth(initialcells = np.array([fp_exact[0],0]))[0]
        g2_exact = g.Growth(initialcells = np.array([1,fp_exact[1]]))[1]/g.Growth(initialcells = np.array([0,fp_exact[1]]))[1]
        g1 = gamma1(1,fp[1],a,y,args.dilution,args.substrate)
        g2 = 1-gamma1(fp[0],1,a,y,args.dilution,args.substrate)
        print '{:6.3f} {:6.3f} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}'.format(a,y,g2*fp[1],g1*fp[0],g2,g1,1-g1_exact,1-g2_exact,fp_exact[0],fp_exact[1])
    print
        
        
        
        
