#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import shapely.geometry as sg

import polygonclasses as pc

def strainparameters(dist,angle):
    return dist*np.cos(angle),dist*np.sin(angle)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1","--infiles_strain1",nargs="*")
    parser.add_argument("-2","--infiles_strain2",nargs="*")
    parser.add_argument("-a","--ExtendToGrowthRates",type=float,default=10)
    parser.add_argument("-W","--WashoutThresholdGrowth",type=float,default=.6)
    parser.add_argument("-Y","--CutAtYield",type=float,default=10)
    parser.add_argument("-s","--step",type=int,default=1)
    parser.add_argument("-D","--baseDilutions",type=float,default=2)
    parser.add_argument("-S","--substrate",type=float,default=1e4)
    
    parser.add_argument("-v","--verbose",action="store_true",default=False)


    parser.add_argument("-d","--dx",type=float,default=.1)
    parser.add_argument("-M","--maxdist",type=float,default=10)
    parser.add_argument("-P","--samplepoints",type=int,default=1000)

    args = parser.parse_args()


    data = pc.Coexistence(args.infiles_strain1,args.infiles_strain2,ExtendToGrowthRates = args.ExtendToGrowthRates, WashoutThresholdGrowth = args.WashoutThresholdGrowth, CutAtYield = args.CutAtYield,verbose = args.verbose)
    baseRegion = data.getPolygon(args.baseDilutions/args.substrate)
    
    alld = np.arange(0,args.maxdist,args.dx,dtype=np.float64)
    
    for d in alld:
        count = 0
        alist,ylist = strainparameters(d,np.random.uniform(2*math.pi,size = args.samplepoints))
        for a,y in zip(alist,ylist):
            p = sg.Point([a,y])
            if baseRegion.contains(p):
                count += 1
        print d,(1.*count)/args.samplepoints
            
            
    
    

if __name__ == "__main__":
    main()















