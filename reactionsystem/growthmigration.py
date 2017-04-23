#!/usr/bin/env python


import numpy as np
import argparse
import sys,math

import reactionsystem as rs

parser = argparse.ArgumentParser()
parser.add_argument("-N","--n0",type=int,default=25)
parser.add_argument("-M","--m0",type=int,default=25)
parser.add_argument("-S","--substrate",type=int,default=1000)
parser.add_argument("-r","--repetitions",type=int,default=10)
parser.add_argument("-m","--mu",type=float,default=1e-3)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-o","--outputsteps",type=int,default=100)
args = parser.parse_args()

r = rs.reactionsystem(indexset="NMRS")

r.add_reaction("N", "M",   rate = args.mu,    coefficients = "N")
r.add_reaction("M", "N",   rate = args.mu,    coefficients = "M")
r.add_reaction("NR", "NN", rate = args.alpha, coefficients = "N")
r.add_reaction("MS", "MM", rate = args.alpha, coefficients = "M")


for rep in range(args.repetitions):
    r.set_population("N",args.n0)
    r.set_population("M",args.m0)
    r.set_population("RS",args.substrate)
    r.set_time(0)
    o=0
    while r.is_present("R"):
        if o%args.outputsteps == 0:
            print "{:8.3f} {:4d} {:4d}".format(r.get_time(),*r.get_populations("NM"))
        o = r.step()
    print "{:8.3f} {:4d} {:4d}".format(r.get_time(),*r.get_populations("NM"))
    print
    
