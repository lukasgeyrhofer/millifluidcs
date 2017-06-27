#!/usr/bin/env python

import numpy as np
import argparse
import sys,math



class inoculumeffect(object):
    def __init__(self,**kwargs):
        self.__generations = np.array(kwargs.get("generations",[10]),dtype=float)
        self.__correlationtime = kwargs.get("correlationtime",10)
        self.__initialpopulation = kwargs.get("initialpopulation",25)
        self.__initialcorrelation = kwargs.get("initialcorrelation",10)
        self.__yieldinterval = np.array([kwargs.get("yieldmin",.5),kwargs.get("yieldmax",1.5)])
        
        
        self.__substrate = 0.5 * (self.__yieldinterval[1] + self.__yieldinterval[0]) * np.power(2.,self.__generations)
        self.__currentrun = 0
        self.__maxrun = len(self.__generations)
        
        self.__coefficient = np.array([np.exp(-1./self.__correlationtime),1. - np.exp(-1./self.__correlationtime)])
        
        self.__seedingpopulation = list()
        x = self.rng()
        for i in range(self.__initialpopulation-1):
            self.__seedingpopulation.append(x)
            for j in range(self.__initialcorrelation):
                x = self.newyield(x)
        
        self.__averagevaluesperbin = 20
        self.__histograms = list()
        self.__finalpopulationsize = list()
        
    def rng(self,count=1):
        return self.__yieldinterval[0] + (self.__yieldinterval[1] - self.__yieldinterval[0]) * np.random.uniform()
            
    def newyield(self,xn):
        return self.__coefficient[0] * xn + self.__coefficient[1] * self.rng()
    
    def run(self):
        if self.__currentrun < self.__maxrun:
            # set initial conditions
            self.__population = list(self.__seedingpopulation[:])
            self.__currentsubstrate = self.__substrate[self.__currentrun]

            # runs until nutrients are out
            while self.add():
                continue
            
            # do statistics on run
            self.__histograms.append(np.histogram(self.__population,range = self.__yieldinterval, bins = int(np.sum(self.__substrate)/(self.__averagevaluesperbin * self.__maxrun))))
            self.__finalpopulationsize.append(len(self.__population))
                                                  

            # update for restart
            self.__currentrun += 1
            self.__seedingpopulation = np.random.choice(self.__population,size = self.__initialpopulation)
            return True
        else:
            return False
    
    def add(self):
        x = self.newyield(np.random.choice(self.__population))
        if self.__currentsubstrate > x:
            self.__currentsubstrate -= x
            self.__population.append(x)
            return True
        else:
            return False
        

    def __getattr__(self,key):
        if key == "finalpopulationsize":
            return self.__finalpopulationsize
        elif key == "histograms":
            if len(self.__histograms) > 0:
                r = np.array([self.__histograms[0][1][:-1] + 0.5 * np.diff(self.__histograms[0][1])])
                for h in self.__histograms:
                    r = np.concatenate([r,np.array([h[0]])],axis=0)
            
                return r.T
            else:
                raise ValueError("no histograms found. run the populations")
        else:
            super(inoculumeffect,self).__getattr__(key)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--generations",default=[10.],type=float,nargs="*")
    parser.add_argument("-t","--correlationtime",default=10,type=float)
    parser.add_argument("-T","--initialcorrelation",type=float,default=10)
    parser.add_argument("-n","--initialpopulation",type=int,default=20)
    parser.add_argument("-y","--yieldmin",type=float,default=.5)
    parser.add_argument("-Y","--yieldmax",type=float,default=1.5)
    
    parser.add_argument("-o","--histo_outfile",default="histogram.txt")
    args = parser.parse_args()


    ie = inoculumeffect(**vars(args))
    
    while ie.run():
        continue
    
    print ie.finalpopulationsize

    np.savetxt(args.histo_outfile,ie.histograms)



if __name__ == "__main__":
    main()




