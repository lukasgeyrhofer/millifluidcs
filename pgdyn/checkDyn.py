#!/usr/bin/env python3


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
        
    
    def dyn(self,time,x):
        #dummy
        return np.zeros(len(x))


    def RungeKutta4(self,xx,tt):
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
        return (" ".join(["{:14.6e}"]*len(self.x))).format(*self.x)


class DynWithPG(Dynamics):
    def __init__(self,**kwargs):
        Dynamics.__init__(self,**kwargs)
        assert len(self.x) == 3
    
    def dyn(self,time,x):
        return np.array( [
            (1. + self.params['eps'] * x[2])*x[0],
            -(1. + self.params['eps'] * x[2])/(1. + self.params['delta'] * x[2]) * x[0],
            self.params['kappa']*x[0]
            ])



class DynDirect(Dynamics):
    def __init__(self,**kwargs):
        Dynamics.__init__(self,**kwargs)
        assert len(self.x) == 2
    
    def dyn(self,time,x):
        return np.array([
            (1+self.params['eps']*x[0])*x[0],
            -(1+self.params['eps']*x[0])/(1+self.params['delta']*x[0])*x[0]
            ])
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e","--epsilon",type=float,default=1e-6)
    parser.add_argument("-d","--delta",type=float,default=1e-6)
    parser.add_argument("-k","--kappa",type=float,default=1e0)
    
    parser.add_argument("-T","--maxtime",type=float,default=20)
    parser.add_argument("-t","--timestep",type=float,default=.1)
    parser.add_argument("-i","--integratetimestep",type=float,default=1e-3)
    
    args = parser.parse_args()
    
    params = {"eps":args.epsilon,"delta":args.delta,"kappa":args.kappa}

    d1 = DynWithPG(step = args.integratetimestep, initialconditions = np.array([1,1e4,0]), params = params)
    d2 = DynDirect(step = args.integratetimestep, initialconditions = np.array([1,1e4]),   params = params)
    
    t = 0
    while t < args.maxtime:
        d1.IntegrateTime(args.timestep)
        d2.IntegrateTime(args.timestep)
        t += args.timestep
        print("{:6.2f}".format(t),d1,d2)
    
    
if __name__ == "__main__":
    main()


    






