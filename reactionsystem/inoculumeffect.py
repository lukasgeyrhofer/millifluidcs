#!/usr/bin/env python

import numpy as np
import argparse
import sys,math



class inoculumeffect(object):
    def __init__(self,**kwargs):
        self.__generations     = kwargs.get("generations",8)
        self.__correlationtime = kwargs.get("correlationtime",8)
        self.__initialpopulation  = kwargs.get("initialpopulation",25)
        self.__initialcorrelation = kwargs.get("initialcorrelation",8)
        self.__initialgenerations = kwargs.get("initialgenerations",8)
        
        self.__yieldinterval = np.array([kwargs.get("yieldmin",.5),kwargs.get("yieldmax",1.5)])
        
        
        #self.__substrate = 0.5 * (self.__yieldinterval[1] + self.__yieldinterval[0]) * np.power(2.,self.__generations) * self.__initialpopulation
        self.__coefficient = np.array([np.exp(-1./self.__correlationtime),1. - np.exp(-1./self.__correlationtime)])
        
        
        self.__histogrambins = 20
        self.__histograms = list()
        self.__finalpopulationsize = list()
        
        self.__haveovernightculture = False
        
    def rng(self,count=1):
        return self.__yieldinterval[0] + (self.__yieldinterval[1] - self.__yieldinterval[0]) * np.random.uniform()
            
    def newyield(self,xn):
        return self.__coefficient[0] * xn + self.__coefficient[1] * self.rng()
    
    
    def run_overnightculture(self,initialpopulation = None, initialcorrelation = None, generations = None):
        if  initialpopulation  is None:
            initialpopulation  = self.__initialpopulation
        if  initialcorrelation is None:
            initialcorrelation = self.__initialcorrelation
        if  generations        is None:
            generations        = self.__initialgenerations
        
        self.__overnightculture = list()
        self.__currentsubstrate = 0.5 * (self.__yieldinterval[1] + self.__yieldinterval[0]) * np.power(2.,generations) * initialpopulation
        
        # make the initial seeding for the overnight cultute
        x = self.rng()
        for i in range(initialpopulation-1):
            self.__overnightculture.append(x)
            # wait 'initialcorrelation' generations before adding a new value, this is only a rough estimate of this distribution
            for j in range(int(initialcorrelation)):
                x = self.newyield(x)
        
        # add more cells, but not to a different population
        while self.add(population = "overnightculture"):
            continue
        
        # we're done here
        self.__haveovernightculture = True
    
    
    def run(self,generations = None,initialpopulation = None):
        # need overnightculture for seeding
        if not self.__haveovernightculture:
            self.run_overnightculture()
            
        # use default values from object creation is no argument given here
        if  initialpopulation is None:
            initialpopulation = self.__initialpopulation
        if  generations       is None:
            generations       = self.__generations
        
        # set initial conditions
        self.__population = list(np.random.choice(self.__overnightculture,size = initialpopulation))
        self.__currentsubstrate = 0.5 * (self.__yieldinterval[1] + self.__yieldinterval[0]) * np.power(2.,generations) * initialpopulation

        # run until nutrients are out
        while self.add():
            continue
        
        # do statistics on run
        self.__histograms.append(np.histogram(self.__population,range = self.__yieldinterval, bins = self.__histogrambins))
        self.__finalpopulationsize.append(len(self.__population))

    
    def add(self,population = "population"):
        x = self.newyield(np.random.choice(self.__dict__["_inoculumeffect__{:s}".format(population)]))
        if self.__currentsubstrate > x:
            self.__currentsubstrate -= x
            self.__dict__["_inoculumeffect__{:s}".format(population)].append(x)
            return True
        else:
            return False
        

    def __getattr__(self,key):
        # reading out those attributes resets them to empty
        if key == "finalpopulationsize":
            fps = self.__finalpopulationsize
            self.__finalpopulationsize = list()
            return fps
        elif key == "histograms":
            if len(self.__histograms) > 0:
                r = np.array([self.__histograms[0][1][:-1] + 0.5 * np.diff(self.__histograms[0][1])])
                for h in self.__histograms:
                    r = np.concatenate([r,np.array([h[0]])],axis=0)
                
                self.__histograms = list()
                return r.T
                
            else:
                raise ValueError("no histograms found. run the populations")
        elif "substraterange":
            return self.__yieldinterval/np.sum(self.__yieldinterval) * np.power(2,self.__generations) * self.__initialpopulation

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--generations",        type=float, default=8.)
    parser.add_argument("-t","--correlationtime",    type=float, default=8.)
    parser.add_argument("-T","--initialcorrelation", type=float, default=8.)
    parser.add_argument("-G","--initialgenerations", type=float, default=8.)
    parser.add_argument("-n","--initialpopulation",  type=int,   default=20)
    
    parser.add_argument("-y","--yieldmin",type=float,default=.5)
    parser.add_argument("-Y","--yieldmax",type=float,default=1.5)
    
    parser.add_argument("-N","--populationcount",type=int,default=1000)
    parser.add_argument("-O","--overnightculturecount",type=int,default=3)
    parser.add_argument("-o","--outfilebasename",default="out")
    args = parser.parse_args()


    substraterange = np.array([args.yieldmin,args.yieldmax])/(args.yieldmin + args.yieldmax) * np.power(2,args.generations) * args.initialpopulation
    ie = inoculumeffect(**vars(args))
    
    
    for i in range(args.overnightculturecount):
        ie.run_overnightculture()
    
        for j in range(args.populationcount):
            ie.run()


        ps,psbin = np.histogram(ie.finalpopulationsize,ie.substraterange,bins = 20)
        psbin = psbin[:-1] + 0.5 * np.diff(psbin)
        np.savetxt("{}-Pdistr{:04d}".format(args.outfilebasename,i),np.transpose(np.array([psbin,ps])))
        np.savetxt("{}-Ydistr{:04d}".format(args.outfilebasename,i),ie.histograms)



if __name__ == "__main__":
    main()




