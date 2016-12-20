#!/usr/bin/env python


import argparse
import numpy as np
import shapely.geometry as sg
import shapely.affinity as sa
import matplotlib.pyplot as plt
import sys
import os.path as op


def rescale(geom,x,y):
    return sa.affine_transform(geom,[x,0,0,y,0,0])

def PlotGraph(axes,polygon,col="#3465a4"):
    if isinstance(polygon,sg.Polygon):
        x,y = polygon.exterior.xy
        axes.plot(x,y,linewidth=2,color=col)
    elif isinstance(polygon,sg.MultiPolygon):
        for p in polygon:
            x,y = p.exterior.xy
            axes.plot(x,y,linewidth=2,color=col)

def PrintGraph(fp,polygon):
    if isinstance(polygon,sg.Polygon):
        xlist,ylist = polygon.exterior.xy
        for x,y in zip(xlist,ylist):
            print >>fp,"{:.6e} {:.6e}".format(x,y)
        print >>fp
    elif isinstance(polygon,sg.MultiPolygon):
        for p in polygon:
            PrintGraph(gp,p)



class Coexistence(object):
    def __init__(self,filenames1,filenames2,ExtendToGrowthRates = 10,CutAtYield = 10,WashoutThresholdGrowth = .6,step = 1,verbose = False):
        
        self.__verbose = verbose
        
        if (filenames1 is None) or (filenames2 is None):
            raise ValueError
        assert len(filenames1) == len(filenames2)
        
        # rawdata
        self.__invasioncurves = np.array([dict(),dict()])

        # store processed data
        self.__polygons    = dict()
        self.__coordinates = dict()
        self.__keys        = list()


        # load data
        for fn in filenames1:
            try:
                self.__invasioncurves[0][self.extractYield(fn)] = np.genfromtxt(fn)
            except:
                pass
        for fn in filenames2:
            try:
                self.__invasioncurves[1][self.extractYield(fn)] = np.genfromtxt(fn)
            except:
                pass
        
        # construct coordinate arrays and polygons
        for key in self.__invasioncurves[0]:
            if not self.__invasioncurves[1].has_key(key):
                raise IndexError
            
            indexCenter1,indexCenter2 = self.getIndexCenter(key)
            direction1 = 1
            direction2 = 1
            if self.__invasioncurves[0][key][indexCenter1 + step,0] < 1 and self.__invasioncurves[0][key][indexCenter1 + step,1] > 1:
                direction1 = -1
            if self.__invasioncurves[1][key][indexCenter2 + step,0] < 1 and self.__invasioncurves[1][key][indexCenter2 + step,1] > 1:
                direction2 = -1
            
            # which is the upper curve?
            slope1 = (self.__invasioncurves[0][key][indexCenter1 + step * direction1,1] - self.__invasioncurves[0][key][indexCenter1,1]) / (self.__invasioncurves[0][key][indexCenter1 + step * direction1,0] - self.__invasioncurves[0][key][indexCenter1,0])
            slope2 = (self.__invasioncurves[1][key][indexCenter2 + step * direction2,1] - self.__invasioncurves[1][key][indexCenter2,1]) / (self.__invasioncurves[1][key][indexCenter2 + step * direction2,0] - self.__invasioncurves[1][key][indexCenter2,0])
            upper = 0
            if slope1 < slope2:
                upper = 1
            
            
            # we know which is the upper curve and have all the raw data
            # construct the polygon out of that
            
            a1 = self.__invasioncurves[upper][key][::-direction1,0]
            y1 = self.__invasioncurves[upper][key][::-direction1,1]
            a2 = self.__invasioncurves[1-upper][key][::direction2,0]
            y2 = self.__invasioncurves[1-upper][key][::direction2,1]
            
            a1 = a1[y1 <= CutAtYield]
            y1 = y1[y1 <= CutAtYield]
            a2 = a2[y2 <= CutAtYield]
            y2 = y2[y2 <= CutAtYield]
            
            a1[a1 <= WashoutThresholdGrowth] = WashoutThresholdGrowth
            a2[a2 <= WashoutThresholdGrowth] = WashoutThresholdGrowth
            
            a = np.concatenate([a1,np.ones(2)*ExtendToGrowthRates,a2,np.array([a2[-1],a1[0]])])
            y = np.concatenate([y1,np.array([y1[-1],y2[0]]),y2,np.ones(2)*y1[0]])
            
            if self.__verbose:
                print >>sys.stderr,"# Load file '{:s}', directions ({:d};{:d}) upper ({:d})".format(key,direction1,direction2,upper)
            self.__coordinates[key] = np.transpose(np.array([a,y]))
            self.__polygons[key] = sg.Polygon(self.__coordinates[key])
            self.__keys.append(float(key))
    
    # small helper routines
    def extractYield(self,filename):
        return (op.basename(filename)).split("_")[1]

    def getIndexCenter(self,key):
        # curves should go through (a,y) = (1,1)
        # find index of closest point to that values
        i1 = ((self.__invasioncurves[0][key][:,0] - 1)**2 + (self.__invasioncurves[0][key][:,1] - 1)**2).argmin()
        i2 = ((self.__invasioncurves[1][key][:,0] - 1)**2 + (self.__invasioncurves[1][key][:,1] - 1)**2).argmin()
        return i1,i2

    def keys(self):
        return list(self.__coordinates.keys())


    # return processed data,
    # if key is a number, get polygon with closest value
    def getPolygon(self,key):
        if isinstance(key,str):
            if self.__polygons.has_key(key):
                return self.__polygons[key]
            else:
                raise KeyError
        elif isinstance(key,(float,np.float,np.float64)):
            mk = np.array([(float(k) - key)**2 for k in self.keys()])
            return self.__polygons[self.keys()[mk.argmin()]]

    # same as above, but with np-coordinate arrays
    def getCoordinates(self,key):
        if isinstance(key,str):
            if self.__coordinates.has_key(key):
                return self.__coordinates[key]
            else:
                raise KeyError
        elif isinstance(key,(float,np.float,np.float64)):
            mk = np.array([(float(k) - key)**2 for k in self.keys()])
            return self.__coordinates[self.keys()[mk.argmin()]]

    # python behavior
    def __int__(self):
        return len(self.keys())
    
    def __str__(self):
        ret = ""
        for key in self.__polygons:
            ret += key + " " + str(self.__polygons[key]) + "\n"
        return ret


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


    data = Coexistence(args.infiles_strain1,args.infiles_strain2,ExtendToGrowthRates = args.ExtendToGrowthRates, WashoutThresholdGrowth = args.WashoutThresholdGrowth, CutAtYield = args.CutAtYield,verbose = args.verbose)

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

        PlotGraph(ax,coexRegion,col="#cc0000")

    PreviousCoexParameters = []
    for sp in StrainParameters:
        if coexRegion.contains(sg.Point([sp[0],sp[1]])):
            curRegion = rescale(data.getPolygon(sp[1]*args.baseDilutions/args.substrate),sp[0],sp[1])
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
                    PlotGraph(ax,curRegion,col="#73d216")
            else:
                if args.verbose:
                    print >> sys.stderr, "# Phasediagram of parameters ( {} : {} ) does not contain reference strain ( 1 : 1 ) in their rescaled coexistence region. Skipping ...".format(*sp)
        else:
            if args.verbose:
                print >> sys.stderr, "# Parameters ( {} : {} ) not contained in coexistence region of reference strain. Cannot compute multistrain phasediagram using these parameters. Skipping ...".format(*sp)

    if args.showGraph:
        PlotGraph(ax,coexRegion)
        plt.show()

    if args.outfile is None:
        fp = sys.stdout
    else:
        try:
            fp = open(args.outfile,"w")
        except:
            fp = sys.stdout
    PrintGraph(fp,coexRegion)
    fp.close()


if __name__ == "__main__":
    main()
