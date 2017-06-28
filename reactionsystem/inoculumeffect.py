#!/usr/bin/env python

import numpy as np
import argparse
import sys,math



class inoculumeffect(object):
    def __init__(self,**kwargs):
        # parse cmdline arguments
        self.__generations          = kwargs.get("generations",8)
        self.__correlation          = kwargs.get("correlation",8)
        self.__seedingsize          = int(kwargs.get("seedingsize",25))
        
        self.__ONinitialcorrelation = kwargs.get("ON_initialcorrelation",8)
        self.__ONgenerations        = kwargs.get("ON_generations",8)
        self.__ONseedingsize        = int(kwargs.get("ON_seedingsize",25))
        
        self.__yieldinterval        = np.array([kwargs.get("yieldmin",.5),kwargs.get("yieldmax",1.5)])
        self.__verbose              = kwargs.get("verbose",False)
        
        # coefficients for faster reference instead of computing them every step
        self.__coefficient          = np.array([np.exp(-1./self.__correlation),1. - np.exp(-1./self.__correlation)])
        
        # statistics, analysis
        self.__histogrambins       = 20
        self.__histograms          = list()
        self.__finalpopulationsize = list()
        
        # have startingconditions?
        self.__haveovernightculture = False
    
    
    
    def rng(self):
        return np.random.uniform(low = self.__yieldinterval[0], high = self.__yieldinterval[1])
            
    def newyield(self,xn):
        return self.__coefficient[0] * xn + self.__coefficient[1] * self.rng()
    
    
    def run_overnightculture(self,seedingsize = None, generations = None, initialcorrelation = None):
        if  seedingsize        is None:
            seedingsize        = self.__ONseedingsize
        if  initialcorrelation is None:
            initialcorrelation = self.__ONinitialcorrelation
        if  generations        is None:
            generations        = self.__ONgenerations
        
        self.__overnightculture = list()
        
        # make the initial seeding for the overnight culture
        x = self.rng()
        for i in range(seedingsize):
            self.__overnightculture.append(x)
            # wait 'initialcorrelation' generations before adding a new value, this is only a rough estimate of this distribution
            for j in range(int(initialcorrelation)):
                x = self.newyield(x)
        
        # from these initial seedings, run on average g generations
        self.__currentsubstrate = np.power(2.,generations) * seedingsize / np.mean(self.__overnightculture)
        # add more cells, but not to a different population
        while self.add(population = "overnightculture"):
            continue
        
        # starting substrate chosen such that the ON culture would take on average g generations to use up all nutrients
        self.__startingsubstrate = np.power(2.,generations) * seedingsize / np.mean(self.__overnightculture)
        # we're done here
        self.__haveovernightculture = True
    
    
    def run(self,seedingsize = None, generations = None):
        # need overnightculture for seeding
        if not self.__haveovernightculture:
            self.run_overnightculture()
            
        # use default values from object creation if no argument given here
        if  seedingsize is None:
            seedingsize = self.__seedingsize
        if  generations is None:
            generations = self.__generations
        
        # set initial conditions
        self.__population = list(np.random.choice(self.__overnightculture,size = seedingsize))
        self.__currentsubstrate = self.__startingsubstrate

        # run until nutrients are out
        while self.add():
            continue
        
        # do statistics on run
        self.__histograms.append(np.histogram(self.__population,range = self.__yieldinterval, bins = self.__histogrambins))
        fps = len(self.__population)
        self.__finalpopulationsize.append(fps)
        return fps

    
    def add(self,population = "population"):
        # use dict representation of self to chose either "self.__population" or "self.__overnightculture"
        x  = self.newyield(np.random.choice(self.__dict__["_inoculumeffect__{:s}".format(population)]))
        xi = 1./x
        if self.__currentsubstrate > xi:
            self.__currentsubstrate -= xi
            self.__dict__["_inoculumeffect__{:s}".format(population)].append(x)
            return True
        else:
            return False
    
    def verbose(self,msg = "", handle = None):
        if self.__verbose:
            if handle is None:
                handle = sys.stdout
            print >>handle, msg
    

    def __getattr__(self,key):
        # reading out those attributes resets them to empty!
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
            return np.sort(np.power(2,self.__generations) * self.__seedingsize / self.__yieldinterval)







def main():
    # run multiple ON cultures with many 'droplets' (=second populations) each
    
    # add cmdline arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--generations",           type = float, default = 8.)
    parser.add_argument("-t","--correlation",           type = float, default = 8.)
    parser.add_argument("-n","--seedingsize",           type = int,   default = 25)
    parser.add_argument("-T","--ON_initialcorrelation", type = float, default = 8.)
    parser.add_argument("-G","--ON_generations",        type = float, default = 8.)
    parser.add_argument("-N","--ON_seedingsize",        type = int,   default = 25)
    
    parser.add_argument("-y","--yieldmin", type = float, default = 0.5)
    parser.add_argument("-Y","--yieldmax", type = float, default = 1.5)
    
    parser.add_argument("-k","--droplets",              type = int, default = 1000)
    parser.add_argument("-O","--overnightculturecount", type = int, default = 3)
    parser.add_argument("-o","--outfilebasename",                   default = "out")
    
    parser.add_argument("-v","--verbose",default=False,action="store_true")
    args = parser.parse_args()

    # initialize object and datastructure
    ie = inoculumeffect(**vars(args))
    
    # loop over different ON cultures
    for i in range(args.overnightculturecount):
        ie.verbose("# starting overnight culture")
        ie.run_overnightculture()
    
        # seed droplets from this ON culture
        for j in range(args.droplets):
            ie.verbose("#   droplet {:4d}".format(j))
            ie.run()

        # reading destroys the data, so only read once
        fps    = ie.finalpopulationsize
        Hyield = ie.histograms
        
        # make histogram for population sizes
        ps,psbin = np.histogram(fps,range = ie.substraterange,bins = 100)
        Hfps = np.transpose(np.array([psbin[:-1] + 0.5 * np.diff(psbin),ps]))
        
        # save histograms to files
        np.savetxt("{}_P{:04d}".format(args.outfilebasename,i),Hfps)
        np.savetxt("{}_Y{:04d}".format(args.outfilebasename,i),Hyield)



if __name__ == "__main__":
    main()




