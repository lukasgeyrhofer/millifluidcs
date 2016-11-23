#!/usr/bin/env python

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


def f(tt,xx):
    return np.array([
        args.growthrate * (1-np.exp(-args.epsilon*x[1])) * x[0],
        x[0]
        ])

parser = argparse.ArgumentParser()
parser.add_argument("-a","--growthrate",type=float,default=1)
parser.add_argument("-e","--epsilon",type=float,default=1e-3)
parser.add_argument("-N","--Ninitial",type=float,default=1)
parser.add_argument("-G","--Ginitial",type=float,default=0)

parser.add_argument("-t","--integrationstep",type=float,default=1e-3)
parser.add_argument("-o","--outputstep",type=float,default=.1)
parser.add_argument("-T","--maxtime",type=float,default=24)
global args
args = parser.parse_args()


x = np.array([args.Ninitial,args.Ginitial])
t = 0

while t < args.maxtime:
    dt = 0
    while dt < args.outputstep:
        x = RungeKutta4(f,x,t+dt,args.integrationstep)
        dt += args.integrationstep
    t += dt
    print "{:14.6f} {:14.6f} {:14.6f}".format(t,*x)
