#!/usr/bin/env python


import numpy as np
import argparse
import sys,math


class Dynamics:
    # General forward integration of dynamics with Runge-Kutta method of 4th order
    # classes have to inherit from here and refine their own 'dyn' function
    def __init__(self,step = 1e-3,requiredpositive = True,initialconditions = None,globaltime = 0,**kwargs):
        self.__step = step
        self.__requiredpositice = requiredpositive

        self.params = kwargs.get('params',None)
        
        if initialconditions is None:
            raise ValueError
        else:
            self.x = np.array(initialconditions)
            
        self.globaltime = globaltime
        
    
    def dyn(time,x):
        #dummy
        return np.zeros(len(x))


    def RungeKutta4(xx,tt):
        # 4th order Runge-Kutta integration scheme
        k1 = self.__step * self.dyn( tt        , xx )
        k2 = self.__step * self.dyn( tt+self.__step/2., xx+k1/2. )
        k3 = self.__step * self.dyn( tt+self.__step/2., xx+k2/2. )
        k4 = self.__step * self.dyn( tt+self.__step   , xx+k3 )
        return xx + (k1+2*k2+2*k3+k4)/6.

    def IntegrateTime(self,time):
        t = 0
        while t <= time:
            self.x = self.RungeKutta4(self.x,self.globaltime + t)
            if self.__requiredpositice:
                self.x[self.x<=0]=0
            t += self.__step
        self.globaltime += time
    
    
    def IntegrateToZero(self,index):
        t = 0
        while self.x[index] > 0:
            self.x = self.RungeKutta4(self.x,self.globaltime + t)
            if self.__requiredpositice:
                self.x[self.x<=0]=0
            t += self.__step
        self.globaltime += t
            

    def __str__(self):
        return (" ".join(["{:10.6f}"]*len(x))).format(*x)


class DynWithPG(Dynamics):
    def __init__(self,**kwargs):
        super(Dynamics,self).__init__(**kwargs)
        assert len(self.x) == 3
    
    def dyn(time,x):
        return np.array( [
            (1. + params['eps'] * x[2])*x[0],
            -(1. + params['eps'] * x[2])/(1. + self.params['delta'] * x[2]) * x[0]
            params['kappa']*x[0]
            ])



class DynDirect(Dynamics):
    def __init__(self,**kwargs):
        super(Dynamics,self).__init__(**kwargs)
        assert len(self.x) == 2
    
    def dyn(time,x):
        return np.array([
            (1+self.params['eps']*x[0])*x[0],
            -(1+self.params['eps']*x[0])/(1+self.params['delta']*x[0])*x[0]
            ])
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e","--epsilon",type=float,default=1e-3)
    parser.add_argument("-d","--delta",type=float,default=1e-3)
    parser.add_argument("-k","--kappa",type=float,default=1e-1)
    
    parser.add_argument("-T","--maxtime",type=float,default=12)
    parser.add_argument("-t","--timestep",type=float,default=.1)
    
    args = parser.parse_args()
    
    params = {"eps":args.epsilon,"delta":args.delta,"kappa":args.kappa}

    d1 = DynWithPG(initialconditions = np.array([1,1e4,0]),params = params)
    d2 = DynDirect(initialconditions = np.array([1,1e4]),params = params)
    
    t = 0
    while t <= args.maxtime:
        d1.IntegrateTime(args.timestep)
        d2.IntegrateTime(args.timestep)
        t += args.timestep
        print t,d1,d2
    
    
if __name__ == "__main__":
    main()


    






