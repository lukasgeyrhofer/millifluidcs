#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ======================================================= #
#  Reaction system to compute growth of microbes          #
#  in a single millifluidic droplet in the presence       #
#  of antibiotics. Both, growth rate and yield factor     #
#  are influenced by increasing antibiotic concentration  #
# ======================================================= #
#  USAGE:                                                 #
#   ./growth.py                                           #
#                                                         #
#  Parameters estimated corresponding to values found in  #
#    Baraban et al., LabChip (2011)                       #
#                                                         #
#  Lukas Geyrhofer, 2016                                  #
# ======================================================= #


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


class dropletdynamics:
    def __init__(self,**kwargs):
        # initialization
        # extract all parameters from the argparser namespace
        self.__start = np.array([kwargs.get('start_bacteria',10.),kwargs.get('start_nutrients',1.5),kwargs.get('start_antibiotics',0)])
        self.__x = self.__start[:]
        self.__antibiotics = {'zmic':kwargs.get('antibiotics_zmic',0.002),'kappa':kwargs.get('antibiotics_kappa',2.),'gamma':kwargs.get('antibiotics_gamma',10)}
        self.__nutrients = {'ks':kwargs.get('nutrients_ks',0.2),'amax':kwargs.get('nutrients_amax',0.02),'yield':kwargs.get('yieldfactor',2e5)}
        self.__time = 0
        self.__steps = 0
        self.__epsilon = kwargs.get('epsilon',0.05)
    
    def __str__(self):
        # return time and all concentrations when object is printed
        s = "{:.3f}".format(self.__time)
        for i in range(len(self.__x)):
            s += " {:e}".format(self.__x[i])
        return s
    
    def dxdt(self,tt,xx):
        # time evolution for system of differential equations
        # it seems that nutrient concentrations and antibiotic concenctration are almost independent effects
        # thus, model growth rate as product of standard Michaelis-Menten kinetics and a sigmoid curve in antibiotic concentration of growth rate (=beta), interpolating between 1 and -gamma (reaching 0 at a concentration of zMIC)
        bzk = np.power(xx[2]/self.__antibiotics['zmic'],self.__antibiotics['kappa'])
        beta = (1-(1+self.__antibiotics['gamma'])*bzk/(bzk + self.__antibiotics['gamma']))
        growthrate  = self.__nutrients['amax'] * xx[1] / (self.__nutrients['ks'] + xx[1]) * beta
        # results in Baraban et al, LabChip(2011) indicate that maximal denisty is also dependent on antibiotic concentration
        # change yield such that growth stops at the same time
        yieldfactor = self.__nutrients['yield'] * beta
        # still unanswered question: why does this depend on initial conditions?
        return np.array([growthrate * xx[0], -growthrate/yieldfactor * xx[0],0])
    
    def step(self):
        # iterate a single step
        self.__x             =  RungeKutta4(self.dxdt,self.__x,self.__time,self.__epsilon)
        # concentrations cannot be smaller than 0
        self.__x[self.__x<0] =  0
        # progress time
        self.__time          += self.__epsilon
        self.__steps         += 1
    

def main():
    parser = argparse.ArgumentParser()
    parser_initialcond = parser.add_argument_group(description = "Initial conditions")
    parser_initialcond.add_argument("-N","--start_bacteria",    type=float, default=10)
    parser_initialcond.add_argument("-S","--start_substrate",   type=float, default=1.5)
    parser_initialcond.add_argument("-B","--start_antibiotics", type=float, default=0)
    
    parser_interaction = parser.add_argument_group(description = "Interaction parameters")
    parser_interaction.add_argument("-z","--antibiotics_zmic",  type=float, default=0.002)
    parser_interaction.add_argument("-k","--antibiotics_kappa", type=float, default=2.)
    parser_interaction.add_argument("-g","--antibiotics_gamma", type=float, default=10.)
    parser_interaction.add_argument("-s","--nutrients_ks",      type=float, default=0.2)
    parser_interaction.add_argument("-a","--nutrients_amax",    type=float, default=0.02)
    parser_interaction.add_argument("-y","--yieldfactor",       type=float, default=2e5)
    
    parser_algorithm   = parser.add_argument_group(description = "Algorithm parameters")
    parser_algorithm.add_argument(  "-e","--epsilon",           type=float, default=0.05)
    parser_algorithm.add_argument(  "-M","--maxsteps",          type=int,   default=2000)
    parser_algorithm.add_argument(  "-O","--outputstep",        type=int,   default=20)
    args = parser.parse_args()
    
    # create droplet object
    droplet = dropletdynamics(**vars(args))
    
    # iterate single steps
    for i in range(args.maxsteps):
        # output ?
        if i%args.outputstep == 0:
            print droplet
        # iterate
        droplet.step()
    
        
if __name__ == "__main__":
    main()
