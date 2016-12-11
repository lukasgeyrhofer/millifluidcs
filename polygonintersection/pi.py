#!/usr/bin/env python


import argparse
import numpy as np
import shapely.geometry as sg
import shapely.affinity as sa
from matplotlib import pyplot
import sys

def rescale(geom,x,y):
    return sa.affine_transform(geom,[1./x,0,0,1./y,0,0])


class Coexistence(object):
    def __init__(self,filenames1,filenames2,ExtendToGrowthRates = 10,CutAtYield = 10,WashoutThresholdGrowth = .6):
        
        if filenames1 is None:
            raise ValueError
        
        assert len(filenames1) == len(filenames2)
        
        self.__invasioncurves = np.array([{},{}])
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
        
        self.__polygons = {}
        self.__keys = list()
        
        for key in self.__invasioncurves[0]:
            if not self.__invasioncurves[1].has_key(key):
                raise IndexError

            indexCenter1,indexCenter2 = self.getIndexCenter(key)
            direction1 = 1
            direction2 = 1
            if self.__invasioncurves[0][key][indexCenter1 + 1,0] < 1:
                direction1 = -1
            if self.__invasioncurves[1][key][indexCenter2 + 1,0] < 1:
                direction2 = -1
            
            # which is the upper curve?
            upper = 0
            if (self.__invasioncurves[0][key][indexCenter1 + direction1,0]-1)/(self.__invasioncurves[0][key][indexCenter1 + direction1,1] - 1) < (self.__invasioncurves[1][key][indexCenter2 + direction2,0] - 1)/(self.__invasioncurves[1][key][indexCenter2 + direction2] - 1):
                upper = 1
            
            
            # we know which is the upper curve and have all the raw data
            # construct the polygon out of that
            
            y1 = self.__invasioncurves[upper][key][::-direction1,0]
            a1 = self.__invasioncurves[upper][key][::-direction1,1]
            y2 = self.__invasioncurves[1-upper][key][::direction2,0]
            a2 = self.__invasioncurves[1-upper][key][::direction2,1]
            
            a1 = a1[y1 <= CutAtYield]
            y1 = y1[y1 <= CutAtYield]
            a2 = a2[y2 <= CutAtYield]
            y2 = y2[y2 <= CutAtYield]
            
            a2[a2<=WashoutThresholdGrowth] = WashoutThresholdGrowth
            
            a = np.concatenate([a1,np.ones(2)*ExtendToGrowthRates,a2])
            y = np.concatenate([y1,np.array([y1[-1],y2[0]]),y2])
            
            self.__polygons[key] = sg.Polygon(np.transpose([a,y]))
            self.__keys.append(float(key))
            
    def extractYield(self,filename):
        return filename.split("_")[1]

    def getIndexCenter(self,key):
        i1 = ((self.__invasioncurves[0][key][:,0] - 1)**2).argmin()
        i2 = ((self.__invasioncurves[0][key][:,0] - 1)**2).argmin()
        return i1,i2



parser = argparse.ArgumentParser()
parser.add_argument("-1","--infiles_strain1",nargs="*")
parser.add_argument("-2","--infiles_strain2",nargs="*")
parser.add_argument("-a","--ExtendToGrowthRates",type=float,default=10)
parser.add_argument("-W","--WashoutThresholdGrowth",type=float,default=.6)
parser.add_argument("-Y","--CutAtYield",type=float,default=10)
args = parser.parse_args()


data = Coexistence(args.infiles_strain1,args.infiles_strain2,ExtendToGrowthRates = args.ExtendToGrowthRates, WashoutThresholdGrowth = args.WashoutThresholdGrowth, CutAtYield = args.CutAtYield)

