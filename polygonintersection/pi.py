#!/usr/bin/env python


import argparse
import numpy as np
import shapely.geometry as sg
import shapely.affinity as sa
from matplotlib import pyplot
import sys

def rescale(geom,x,y):
    return sa.affine_transform(geom,[1./x,0,0,1./y,0,0])


def 


def loadBoundaries(file1,file2):
    try:
        data1 = np.genfromtxt(file1)
        data2 = np.genfromtxt(file2)
    except:
        return None
    y1 = data1[:,0]
    a1 = data1[:,1]
    y2 = data2[:,0]
    a2 = data2[:,1]
    index1_1 = (y1-1).argmin()
    index1_2 = (y2-1).argmin()
    direction_curve1 = 1 y1[
    return sg.Polygon()




parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs=2)
parser.add_argument("-a","--aratio_washout",type=float,default=.1)
args = parser.parse_args()


try:
    data1 = np.genfromtxt(args.infiles[0])
    data2 = np.genfromtxt(args.infiles[1])
except:
    print >> sys.stderr,"could not open files"
    exit(1)


focalStrain = sg.Point([1,1])


a = loadBoundaries("as","as")
exit(1)

branch1u = data1[data1[:,0] > 1]

fig = pyplot.figure()
ax = fig.add_subplot(111)


t = 2.*np.pi*np.arange(10)/10.
coords = np.array([np.sin(t),np.cos(t)]).transpose()
coordsshift = coords + np.array([.5,.5])

polygon1 = sg.Polygon(coords)
polygon2 = sg.Polygon(coordsshift)


plot_line(ax,polygon1.exterior)
plot_line(ax,polygon2.exterior)
#plot_line(ax,polygon1.exterior)

plot_line2(ax,polygon1.intersection(polygon2).exterior)

two = sg.Point([2,2])

tmp = polygon1.intersection(two)

print polygon1.contains(two)


#pyplot.show()
