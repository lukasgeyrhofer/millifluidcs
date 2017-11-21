#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

sys.path.append(sys.path[0] + '/../')

import growthclasses as gc
import pickle


def interpolate(mat1,mat2,weight1 = .5,weight2 = .5):
    if weight1 >= 0 and weight2 >= 0 and weight1 + weight2 > 0:
        return weight1/(weight1 + weight2) * mat1 + weight2/(weight1 + weight2) * mat2
    else:
        return None

def getWeights(a,b,c):
    if a - b != 0:
        return (a - c)/(a - b), (c-b)/(a-b)
    
def deviation(mat1,mat2):
    return np.sqrt(np.sum((mat1 - mat2)*(mat1 - mat2)))

def deviation_relative(mat1,mat2):
    return np.sqrt(np.sum(np.nan_to_num((mat1 - mat2)*(mat1 - mat2)/mat1/mat1)))

parser = argparse.ArgumentParser()
parser.add_argument("-A","--gmA")
parser.add_argument("-a","--valueA",type=float,default=1)
parser.add_argument("-B","--gmB")
parser.add_argument("-b","--valueB",type=float,default=2)
parser.add_argument("-C","--gmC")
parser.add_argument("-c","--valueC",type=float,default=1.5)
args = parser.parse_args()

try:
    gA = pickle.load(open(args.gmA))
    gB = pickle.load(open(args.gmB))
    gC = pickle.load(open(args.gmC))
except:
    raise IOError("could not load input files ('{}', '{}', '{}')".format(args.gmA,args.gmB,args.gmC))

#print gA.growthmatrixgrid

assert np.any(gA.growthmatrixgrid[0] == gB.growthmatrixgrid[0])
assert np.any(gA.growthmatrixgrid[1] == gB.growthmatrixgrid[1])
assert np.any(gA.growthmatrixgrid[0] == gC.growthmatrixgrid[0])
assert np.any(gA.growthmatrixgrid[1] == gC.growthmatrixgrid[1])

a,b = getWeights(args.valueA,args.valueB,args.valueC)

intpC0 = interpolate(gA.growthmatrix[:,:,0],gB.growthmatrix[:,:,0],a,b)
intpC1 = interpolate(gA.growthmatrix[:,:,1],gB.growthmatrix[:,:,1],a,b)

dev0 = deviation_relative(gC.growthmatrix[:,:,0],intpC0)
dev1 = deviation_relative(gC.growthmatrix[:,:,1],intpC1)

print "{:.3e} {:.3e}".format(dev0,dev1)








