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
            raise ValueError,"filenames empty!"
        assert len(filenames1) == len(filenames2)
        
        # rawdata
        self.__invasioncurves = np.array([dict(),dict()])

        # store processed data
        self.__polygons    = dict()
        self.__keys        = list()


        self.ExtendToGrowthRates = ExtendToGrowthRates
        self.WashoutThresholdGrowth = WashoutThresholdGrowth
        self.CutAtYield = CutAtYield

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
            
            
            if self.__verbose:
                print >>sys.stderr,"# Load file '{:s}', directions ({:d};{:d}) upper ({:d})".format(key,direction1,direction2,upper)
            self.__polygons[key] = self.makePolygon(self.__invasioncurves[0][key],self.__invasioncurves[1][key])
            self.__keys.append(float(key))
        
        self.__keys = np.sort(self.__keys)
    
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
        return self.__keys


    def makePolygon(self,contour1,contour2):
        indexCenter1 = ((contour1[:,0] - 1)**2 + (contour1[:,1] - 1)**2).argmin()
        indexCenter2 = ((contour2[:,0] - 1)**2 + (contour2[:,1] - 1)**2).argmin()

        direction1 = 1
        direction2 = 1
        if contour1[indexCenter1 + step,0] < 1 and contour1[indexCenter1 + step,1] > 1:
            direction1 = -1
        if contour2[indexCenter2 + step,0] < 1 and contour2[indexCenter2 + step,1] > 1:
            direction2 = -1
        
        # which is the upper curve?
        slope1 = (contour1[indexCenter1 + step * direction1,1] - contour1[indexCenter1,1]) / (contour1[indexCenter1 + step * direction1,0] - contour1[indexCenter1,0])
        slope2 = (contour2[indexCenter2 + step * direction2,1] - contour2[indexCenter2,1]) / (contour2[indexCenter2 + step * direction2,0] - contour2[indexCenter2,0])
        
        if slope1 < slope2:
            a1 = contour2[::-direction1,0]
            y1 = contour2[::-direction1,1]
            a2 = contour1[::direction2,0]
            y2 = contour1[::direction2,1]
        else:
            a1 = contour2[::-direction1,0]
            y1 = contour2[::-direction1,1]
            a2 = contour1[::direction2,0]
            y2 = contour1[::direction2,1]
        
        
        # we know which is the upper curve and have all the raw data
        # construct the polygon out of that
        
        
        a1 = a1[y1 <= self.CutAtYield]
        y1 = y1[y1 <= self.CutAtYield]
        a2 = a2[y2 <= self.CutAtYield]
        y2 = y2[y2 <= self.CutAtYield]
        
        a1[a1 <= self.WashoutThresholdGrowth] = self.WashoutThresholdGrowth
        a2[a2 <= self.WashoutThresholdGrowth] = self.WashoutThresholdGrowth
        
        a = np.concatenate([a1,np.ones(2)*self.ExtendToGrowthRates,a2,np.array([a2[-1],a1[0]])])
        y = np.concatenate([y1,np.array([y1[-1],y2[0]]),y2,np.ones(2)*y1[0]])
        return sg.Polygon(np.transpose(np.array([a,y])))

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
    def getCoordinates(self,key,multipolygonindex = 0):
        return np.array(sg.mapping(self.getPolygon(key))['coordinates'])

    # python behavior
    def __int__(self):
        return len(self.keys())
    
    def __str__(self):
        ret = ""
        for key in self.__polygons:
            ret += key + " " + str(self.__polygons[key]) + "\n"
        return ret

