#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc


def iterationmap(r,p,gm1,gm2):
    n1 = (1-p)*r*gm1[1,0] + p*r*r*gm1[2,0] + p*r*(1-r)*gm1[1,1]
    n2 = (1-p)*(1-r)*gm2[0,1] + p*(1-r)*(1-r)*gm2[0,2] + p*r*(1-r)*gm2[1,1]
    return n1/(n1+n2)

def diterationmap(r,p,gm1,gm2):
    n1  = (1-p)*r*gm1[1,0] + p*r*r*gm1[2,0] + p*r*(1-r)*gm1[1,1]
    dn1 = (1-p)*gm1[1,0] + 2*p*r*gm1[2,0] + p*(1-2*r)*gm1[1,1]
    n2  = (1-p)*(1-r)*gm2[0,1] + p*(1-r)*(1-r)*gm2[0,2] + p*r*(1-r)*gm2[1,1]
    dn2 = -(1-p)*gm2[0,1] - 2*(1-r)*p*gm2[0,2] + p*(1-2*r)*gm2[1,1]
    return dn1/(n1+n2) + n1*(dn1+dn2)/(n1+n2)**2

parser = argparse.ArgumentParser()
parser = gc.AddGrowthParameters(parser)
parser.add_argument("-p","--imperfectseeding",type=float,default=.2)
parser.add_argument("-D","--delta",type=float,default=.2)
parser.add_argument("-N","--NRprecison",type=float,default=1e-10)
parser.add_argument("-A","--NRalpha",type=float,default=1.)
args = parser.parse_args()


#print growth1,growth2

for yexp in np.arange(-1,1,args.delta):
    for aexp in np.arange(-1,1,args.delta):
        y = np.array([10**yexp,1.])
        a = np.array([10**aexp,1.])
        #print y,a
        g = gc.GrowthDynamics(growthrates = a,yieldrates = y,dilution = 1, mixingtime = args.mixingtime, substrate = args.substrateconcentration)
        growth1,growth2 = g.getGrowthMatrix(3)
        
        curr = 0.5
        lastr = 1
        while not (-args.NRprecison < curr - lastr < args.NRprecison):
            lastr = curr
            #curr -= args.NRalpha * iterationmap(lastr,args.imperfectseeding,growth1,growth2) / diterationmap(lastr,args.imperfectseeding,growth1,growth2)
            curr = iterationmap(lastr,args.imperfectseeding,growth1,growth2)
        
        print "{:.4f} {:.4f} {:.4f} {:.6f} {:9.6f}".format(y[0],a[0],args.imperfectseeding,curr,diterationmap(curr,args.imperfectseeding,growth1,growth2))
    print
        
        
        
