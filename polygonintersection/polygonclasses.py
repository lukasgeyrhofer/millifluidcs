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
    def __init__(self,filenames1,filenames2,ExtendToGrowthRates = 10,CutAtYield = 10,WashoutThresholdGrowth = .6,step = 1, DirectReturnThreshold = 1e-100, verbose = False):
        
        self.__verbose = verbose
        
        if (filenames1 is None) or (filenames2 is None):
            raise ValueError,"filenames empty!"
        assert len(filenames1) == len(filenames2), "each polygon needs two countour files"
        
        # rawdata
        self.__invasioncurves = np.array([dict(),dict()])

        # store processed data
        self.__polygons    = list()
        self.__keys        = list()
        self.__indices     = dict()


        self.ExtendToGrowthRates = ExtendToGrowthRates
        self.WashoutThresholdGrowth = WashoutThresholdGrowth
        self.CutAtYield = CutAtYield
        self.step = step
        self.DirectReturnThreshold = DirectReturnThreshold

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
                raise IndexError,"Countour keys do not match"
            
            
            self.__polygons.append( self.makePolygon(self.__invasioncurves[0][key],self.__invasioncurves[1][key]) )
            self.__keys.append(key)
            self.__indices[key] = len(self.__keys)-1
            if self.__verbose:
                print >>sys.stderr,"# Generated polygon with key '{:s}'".format(key)
        
        self.__keysfloat = np.array([float(k) for k in self.__keys])
        self.__keysindex = np.transpose(np.array([np.arange(len(self.__keys)),self.__keysfloat,np.log(self.__keysfloat) ] ))
        
        #print self.__keysindex
    
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
        if contour1[indexCenter1 + self.step,0] < 1 and contour1[indexCenter1 + self.step,1] > 1:
            direction1 = -1
        if contour2[indexCenter2 + self.step,0] < 1 and contour2[indexCenter2 + self.step,1] > 1:
            direction2 = -1
        
        # which is the upper curve?
        slope1 = (contour1[indexCenter1 + self.step * direction1,1] - contour1[indexCenter1,1]) / (contour1[indexCenter1 + self.step * direction1,0] - contour1[indexCenter1,0])
        slope2 = (contour2[indexCenter2 + self.step * direction2,1] - contour2[indexCenter2,1]) / (contour2[indexCenter2 + self.step * direction2,0] - contour2[indexCenter2,0])
        
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


    def dist(self,Point1,ListPoints2):
        # log distance in fitness space
        return np.array([np.linalg.norm(np.log(Point1)-np.log(Point2)) for Point2 in ListPoints2])

    def interpolateContour(self,contour1,contour2,mixing = .5):
        # start with first entry
        newcontour = np.array([mixing * contour1[0] + (1-mixing) * contour2[self.dist(contour1[0],contour2).argmin()]])
        # then iterate over others
        for Point1 in contour1[1:]:
            newcontour = np.concatenate([newcontour,np.array([mixing * Point1 + (1-mixing) * contour2[self.dist(Point1,contour2).argmin()]])])
        return newcontour

    # return processed data,
    # if key is a number, get polygon with closest value
    def getPolygon(self,key):
        if isinstance(key,str):
            if key in self.__keys:
                return self.__polygons[np.where(self.__keys == key)]
            else:
                raise KeyError
        elif isinstance(key,(float,np.float,np.float64)):
            # check first before starting expensive calculations
                
            indexdata = np.transpose(np.array([self.__keysindex[:,0],(self.__keysindex[:,1] - key)**2]) )
            indexdata = indexdata[indexdata[:,1].argsort()]
            
            closestIndex = int(indexdata[0,0])
            secondIndex  = int(indexdata[1,0])
            
            if indexdata[0,1] < self.DirectReturnThreshold:
                return self.__polygons[closestIndex]
            else:
                #print closestIndex,secondIndex
            
                mixing = (np.log(key) - self.__keysindex[secondIndex,2]) / (self.__keysindex[closestIndex,2] - self.__keysindex[secondIndex,2])
                contour1 = self.interpolateContour(self.__invasioncurves[0][self.__keys[closestIndex]],self.__invasioncurves[0][self.__keys[secondIndex]],mixing)
                contour2 = self.interpolateContour(self.__invasioncurves[1][self.__keys[closestIndex]],self.__invasioncurves[1][self.__keys[secondIndex]],mixing)
                
                return self.makePolygon(contour1,contour2)
            
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

