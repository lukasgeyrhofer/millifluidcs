#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import shapely.geometry as sg
import shapely.affinity as sa

import matplotlib.pyplot  as plt
import matplotlib.patches as patches

import polygonclasses as pc



def samplepoints(count = 10000,MaxVal = np.array([2,2])):
    if MaxVal[0] > 1 and MaxVal[1] > 1:
        return np.array([np.exp(np.random.uniform(0,np.log(MaxVal[0]),size=count)),np.exp(np.random.uniform(0,np.log(MaxVal[1]),size=count))]).transpose()
    else:
        return None



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

    parser.add_argument("-P","--StrainParameters",type=float,nargs="*",default=[1.2,0.65])
    parser.add_argument("-M","--maxsamples",type=int,default=10000)
    args = parser.parse_args()
    
    data = pc.Coexistence(args.infiles_strain1,args.infiles_strain2,ExtendToGrowthRates = args.ExtendToGrowthRates, WashoutThresholdGrowth = args.WashoutThresholdGrowth, CutAtYield = args.CutAtYield,verbose = args.verbose)
    
    baseRegion = data.getPolygon(args.baseDilutions/args.substrate)
    minA,minY,maxA,maxY = baseRegion.bounds
    
    StrainParameters = []
    i = 0
    if not args.StrainParameters is None:
        while i < len(args.StrainParameters):
            try:
                if minA < args.StrainParameters[i] < maxA and minY < args.StrainParameters[i+1] < maxY:
                    StrainParameters.append(np.array([args.StrainParameters[i],args.StrainParameters[i+1]]))
            except:
                pass
            i += 2

    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.set_yscale("log")


    pc.PlotGraph(ax,baseRegion)

    for ProdStrainParam in StrainParameters:
        prodStrain = sg.Point(ProdStrainParam)
        inside = np.array([ProdStrainParam])
        outside = np.array([ProdStrainParam])
        if baseRegion.contains(prodStrain):
            if ProdStrainParam[0] < 1 and ProdStrainParam[1] > 1:
                MaxVal = np.array([ProdStrainParam[0]/minA,ProdStrainParam[1]])
            else:
                MaxVal = np.array([ProdStrainParam[0],ProdStrainParam[1]/minY])
            for sample in samplepoints(args.maxsamples,MaxVal):
                
                newRegion = data.getPolygon(args.baseDilutions/args.substrate * sample[1])
                
                ProdStrain = sg.Point(ProdStrainParam/sample)
                if newRegion.contains(ProdStrain):
                    if args.verbose:
                        print "inside",sample
                    inside = np.concatenate([inside,np.array([sample])])
                else:
                    if args.verbose:
                        print "outside",sample
                    outside = np.concatenate([outside,np.array([sample])])
            
            ax.add_patch(patches.Rectangle((1,1),MaxVal[0]-1,MaxVal[1]-1,facecolor='None'))
            ax.scatter(inside[:,0],inside[:,1],s=3,zorder=2,c='Green')
                         
                    
            
    print len(inside),len(outside),args.maxsamples
    plt.show()
    


if __name__ == "__main__":
    main()
