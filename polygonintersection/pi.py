#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import polygonclasses as pc

def RandomSamplePoints(polygon,inside = True,count = 1000):
    samplepoints = list()
    x0,y0,x1,y1 = polygon.bounds
    while len(samplepoints) < count:
        x = np.random.uniform(x0,x1)
        y = np.random.uniform(y0,y1)
        p = sg.Point([x,y])
        if inside:
            if polygon.contains(p):
                samplepoints.append(np.array([x,y]))
        else:
            samplepoints.append(np.array([x,y]))
    return samplepoints

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
    
    parser.add_argument("-P","--StrainParameters",type=float,nargs="*",default=None)
    parser.add_argument("-o","--outfile",default=None)
    parser.add_argument("-G","--showGraph",default=False,action="store_true")
    args = parser.parse_args()


    data = pc.Coexistence(args.infiles_strain1,args.infiles_strain2,ExtendToGrowthRates = args.ExtendToGrowthRates, WashoutThresholdGrowth = args.WashoutThresholdGrowth, CutAtYield = args.CutAtYield,verbose = args.verbose)

    coexRegion = data.getPolygon(args.baseDilutions/args.substrate)
    baseRegion = data.getPolygon(args.baseDilutions/args.substrate)

    centerPoint = sg.Point([1,1])

    miA,miY,maA,maY = coexRegion.bounds

    StrainParameters = []
    i = 0
    if not args.StrainParameters is None:
        while i < len(args.StrainParameters):
            try:
                StrainParameters.append(np.array([args.StrainParameters[i],args.StrainParameters[i+1]]))
            except:
                pass
            i += 2

    if args.showGraph:
        fig = plt.figure(1, figsize=(5,5), dpi=90)
        ax = fig.add_subplot(111)
        ax.set_xscale("log")
        ax.set_yscale("log")

        pc.PlotGraph(ax,coexRegion,col="#cc0000")

    PreviousCoexParameters = []
    for sp in StrainParameters:
        if coexRegion.contains(sg.Point([sp[0],sp[1]])):
            curRegion = pc.rescale(data.getPolygon(sp[1]*args.baseDilutions/args.substrate),sp[0],sp[1])
            if curRegion.contains(centerPoint):
                containsOthers = True
                for pp in PreviousCoexParameters:
                    if not curRegion.contains(pp):
                        containsOthers = False
                if containsOthers:
                    coexRegion = coexRegion.intersection(curRegion)
                    PreviousCoexParameters.append(sg.Point([sp[0],sp[1]]))
                else:
                    if args.verbose:
                        print >> sys.stderr, "# Phasediagram of parameters ( {} : {} ) does not intersect with all previous strains. Skipping ...".format(*sp)
                    if args.showGraph:
                        pc.PlotGraph(ax,curRegion,col = "#d3d7cf")
                if args.showGraph:
                    pc.PlotGraph(ax,curRegion,col="#73d216")
            else:
                if args.verbose:
                    print >> sys.stderr, "# Phasediagram of parameters ( {} : {} ) does not contain reference strain ( 1 : 1 ) in their rescaled coexistence region. Skipping ...".format(*sp)
                if args.showGraph:
                    pc.PlotGraph(ax,curRegion,col = "#d3d7cf")
        else:
            if args.verbose:
                print >> sys.stderr, "# Parameters ( {} : {} ) not contained in coexistence region of reference strain. Cannot compute multistrain phasediagram using these parameters. Skipping ...".format(*sp)
            if args.showGraph:
                pc.PlotGraph(ax,curRegion,col = "#d3d7cf")

    if args.showGraph:
        pc.PlotGraph(ax,coexRegion)
        plt.show()

    if args.outfile is None:
        fp = sys.stdout
    else:
        try:
            fp = open(args.outfile,"w")
        except:
            fp = sys.stdout
    pc.PrintGraph(fp,coexRegion)
    fp.close()


if __name__ == "__main__":
    main()
