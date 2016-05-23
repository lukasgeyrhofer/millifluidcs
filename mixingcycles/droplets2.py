#!/usr/bin/env python


import numpy as np
import argparse
import sys,math
from scipy import stats

from growthclasses import growthdynamics



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n","--droplets",type=int,default=100)
    parser.add_argument("-d","--dilution",type=float,default=1e-4)
    parser.add_argument("-S","--substrate",type=float,default=1e5)
    parser.add_argument("-g","--growthrates",type=float,nargs="*",default=[1,2])
    parser.add_argument("-y","--yieldrates",type=float,nargs="*",default=[2,1])
    parser.add_argument("-N","--initialpopulation",type=float,nargs="*",default=[1e6,1e6])
    parser.add_argument("-T","--mixingtime",type=float,default=1e5)
    parser.add_argument("-M","--maxsteps",type=int,default=10000)
    parser.add_argument("-O","--outputsteps",type=int,default=100)
    parser.add_argument("-v","--verbose",action="store_true",default=False)
    parser.add_argument("-q","--debug",action="store_true",default=False)

    args = parser.parse_args()
    
    g = growthdynamics(growthrates = np.array(args.growthrates), yieldrates = np.array(args.yieldrates), mixingtime = args.mixingtime, dilution = args.dilution, substrate = args.substrate)
    pool = np.array(args.initialpopulation)
    #pool *= args.dilution/args.droplets
    
    for m in range(args.maxsteps):
        
        popdrop = np.random.poisson(pool,size = (args.droplets,g.numstrains))
        seeding = np.mean(popdrop,axis=0)
        pool = np.zeros(g.numstrains)
        
        if m%args.outputsteps==0:
            print "{:5d} {:8.3f} {:8.3f}".format(m,seeding[0],seeding[1])
        
        for k in range(args.droplets):
            pool += g.getGrowth(np.array(popdrop[k]))/(1.*args.droplets)



if __name__ == "__main__":
    main()
