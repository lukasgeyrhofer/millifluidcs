#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math



def RungeKutta4(func,xx,tt,step):
  # 4th order Runge-Kutta integration scheme
  k1 = step * func( tt        , xx )
  k2 = step * func( tt+step/2., xx+k1/2. )
  k3 = step * func( tt+step/2., xx+k2/2. )
  k4 = step * func( tt+step   , xx+k3 )
  return xx + (k1+2*k2+2*k3+k4)/6.


class growthdynamics:
    def __init__(self,growthtype="monod",**kwargs):
        self.__nutrients = {'ks':kwargs.get('halfgrowthfactor'), 'amax':kwargs.get('maxgrowthrate'), 'yield':kwargs.get('yieldfactor')}
        self.__popsizeparameters = {'logmin':np.log10(kwargs.get('popsizemin')),'logmax':np.log10(kwargs.get('popsizemax')),'logstep':kwargs.get('stepfactorpopsize')}
        self.__substrateparameters = {'logmin':np.log10(kwargs.get('substratemin')),'logmax':np.log10(kwargs.get('substratemax')),'logstep':kwargs.get('stepfactorsubstrate')}
        
        self.__popsize   = np.power(10.,np.arange(self.__popsizeparameters['logmin'],self.__popsizeparameters['logmax']+self.__popsizeparameters['logstep'],self.__popsizeparameters['logstep']))
        self.__substrate = np.power(10.,np.arange(self.__substrateparameters['logmin'],self.__substrateparameters['logmax']+self.__substrateparameters['logstep'],self.__substrateparameters['logstep']))
        self.__stoptimes = np.full(self.shape,None)
        
        self.__epsilon = kwargs.get("epsilon")
        
        if growthtype == "monod":
            self.__growthfunction = self.growthMonod
        elif growthtype == "max":
            self.__growthfunction = self.growthMax
        elif growthtype == "linear":
            self.__growthfunction = self.growthLinear
        else:
            raise IndexError
        
        
    def __getattr__(self,key):
        if key == "popsize":
            return self.__popsize
        elif key == "substrate":
            return self.__substrate
        elif key == "shape":
            return (len(self.__popsize),len(self.__substrate))
        elif key == "stoptimes":
            return self.__stoptimes
        else:
            raise KeyError

    def growthMonod(self,t,x):
        growthrate = self.__nutrients['amax'] * x[1]/(self.__nutrients['ks'] + x[1])
        return np.array([growthrate * x[0], - growthrate/self.__nutrients['yield'] * x[0]])
    
    def growthMax(self,t,x):
        return np.array([self.__nutrients['amax'] * x[0], - self.__nutrients['amax']/self.__nutrients['yield'] * x[0]])
    
    def growthLinear(self,t,x):
        return np.array([self.__nutrients['amax']/self.__nutrients['ks'] * x[0]*x[1] , -self.__nutrients['amax']/(self.__nutrients['yield']*self.__nutrients['ks']) * x[0] * x[1]])


    def computeStopTimes(self):
        for i in range(len(self.__popsize)):
            for j in range(len(self.__substrate)):
                x = np.array([self.popsize[i],self.substrate[j]])
                t = 0
                while x[1] > 1./self.__nutrients['yield']:
                    x = RungeKutta4(self.__growthfunction,x,t,self.__epsilon)
                    t += self.__epsilon
                self.__stoptimes[i,j] = t
            





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n","--popsizemin",type=float,default=1)
    parser.add_argument("-N","--popsizemax",type=float,default=1000)
    parser.add_argument("-s","--substratemin",type=float,default=0.1)
    parser.add_argument("-S","--substratemax",type=float,default=10)
    parser.add_argument("-i","--stepfactorpopsize",type=float,default=0.5)
    parser.add_argument("-I","--stepfactorsubstrate",type=float,default=0.5)
    parser.add_argument("-y","--yieldfactor",type=float,default=2e5)
    parser.add_argument("-a","--maxgrowthrate",type=float,default=0.02)
    parser.add_argument("-k","--halfgrowthfactor",type=float,default=0.2)
    parser.add_argument("-e","--epsilon",type=float,default=0.01)
    args = parser.parse_args()

    growthmonod = growthdynamics("monod", **vars(args))
    growthmonod.computeStopTimes()
    growthmax   = growthdynamics("max",   **vars(args))
    growthmax.computeStopTimes()
    growthlin   = growthdynamics("linear",**vars(args))
    growthlin.computeStopTimes()

    for i in range(growthmonod.shape[0]):
        for j in range(growthmonod.shape[1]):
            print("{:.4e} {:.4e} {:.4e} {:.4e} {:.4e}".format(growthmonod.popsize[i],growthmonod.substrate[j],growthmonod.stoptimes[i,j],growthmax.stoptimes[i,j],growthlin.stoptimes[i,j]))
        print
    
    

if __name__ == "__main__":
    main()








