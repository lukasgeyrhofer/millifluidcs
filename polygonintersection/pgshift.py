#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import shapely.geometry as sg
import shapely.affinity as sa

import polygonclasses as pc


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

    parser.add_argument("-P","--strainparameters",type=float,nargs="*",default=[1.2,0.8])
    args = parser.parse_args()
    
    data = pc.Coexistence(args.infiles_strain1,args.infiles_strain2,ExtendToGrowthRates = args.ExtendToGrowthRates, WashoutThresholdGrowth = args.WashoutThresholdGrowth, CutAtYield = args.CutAtYield,verbose = args.verbose)
    
    
    baseRegion = data.getPolygon(args.baseDilutions/args.substrate)
    
    StrainParameters = []
    i = 0
    if not args.StrainParameters is None:
        while i < len(args.StrainParameters):
            try:
                StrainParameters.append(np.array([args.StrainParameters[i],args.StrainParameters[i+1]]))
            except:
                pass
            i += 2


    for s in StrainParameters:
        prodStrain = sg.Point(s)
        if baseRegion.contains(prodStrain):
            print "yey",s




if __name__ == "__main__":
    main()
