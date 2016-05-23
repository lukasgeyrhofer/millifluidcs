#!/usr/bin/env python


import numpy as np
import argparse
import sys,math
from scipy import stats

from growthclasses import growthdynamics



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-k","--droplets",type=int,default=100)
    parser.add_argument("-a","--growthrates",nargs=2,default=np.array([1.,1.]))
    parser.add_argument("-y","--yieldrates",nargs=2,default=np.array([.90,1.]))
    parser.add_argument("-s","--substrate",default=1e5,type=float)
    parser.add_argument("-d","--dilution",type=float,default=2e-4)
    parser.add_argument("-T","--mixingtime",type=float,default=24.)

    parser.add_argument("-N","--initialsize",nargs=2,default=np.array([1e5,1e5]))
    parser.add_argument("-m","--mixingcycles",type=int,default=20)


    parser.add_argument("-O","--outputsteps",type=int,default=100)
    parser.add_argument("-v","--verbose",action="store_true",default=False)
    parser.add_argument("-q","--debug",action="store_true",default=False)

    args = parser.parse_args()
    
    g = growthdynamics(growthrates = np.array(args.growthrates,dtype = float), yieldrates = np.array(args.yieldrates,dtype = float), mixingtime = args.mixingtime, dilution = args.dilution, substrate = args.substrate)
    pool = np.array(args.initialsize,dtype = float) * g.dilution / args.droplets
    
    for m in range(args.mixingcycles):
        
        popdrop = np.random.poisson(pool,size = (args.droplets,g.numstrains))
        seeding = np.mean(popdrop,axis=0)
        pool = np.zeros(g.numstrains)
        
        if m%args.outputsteps==0:
            print "{:5d} {:8.3f} {:8.3f}".format(m,seeding[0],seeding[1])
        
        for k in range(args.droplets):
            pool += g.getGrowth(np.array(popdrop[k]))/(1.*args.droplets)



if __name__ == "__main__":
    main()
