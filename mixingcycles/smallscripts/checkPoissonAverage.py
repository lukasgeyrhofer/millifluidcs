#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + "/..")
import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-M","--maxM",type=int,default=300)

parser_seed = parser.add_argument_group(description = "Average Seeding")
parser_seed.add_argument("-n","--minN",type=float,default=.1)
parser_seed.add_argument("-N","--maxN",type=float,default=10)
parser_seed.add_argument("-l","--logscaleN",default=False,action="store_true")
parser_seed.add_argument("-k","--stepsN",type=int,default=91)

parser_exp = parser.add_argument_group(description = "Exponents")
parser_exp.add_argument("-a","--minA",type=float,default=.1)
parser_exp.add_argument("-A","--maxA",type=float,default=2)
parser_exp.add_argument("-K","--stepsA",type=int,default=20)
parser_exp.add_argument("-L","--logscaleA",default=False,action="store_true")
parser_exp.add_argument("-m","--mirrorA",default=False,action="store_true")

args = parser.parse_args()


if args.logscaleA:
    alist = np.exp(np.linspace(start = np.log(args.minA), stop = np.log(args.maxA), num = args.stepsA))
else:
    alist = np.linspace(start = args.minA, stop = args.maxA, num = args.stepsA)
if args.mirrorA:
    alist = np.concatenate([-alist[::-1],alist])
if args.logscaleN:
    nlist = np.exp(np.linspace(start = np.log(args.minN), stop = np.log(args.maxN), num = args.stepsN))
else:
    nlist = np.linspace(start = args.minN, stop = args.maxN, num = args.stepsN)

m = np.arange(args.maxM,dtype=float)

mlogm = np.arange(args.maxM,dtype=float)
mlogm[m>0] *= np.log(m[m>0])

for n in nlist:
    px = gc.PoissonSeedingVectors(m,[n])
    for a in alist:
        ma = np.power(m[m>0],a)

        ma = np.nan_to_num(ma)
        
        res1   = np.dot(ma,px[0][m>0])
        appr1  = np.power(n,a)
        appr2  = 1 + a * np.log(n)

        #eps1 = np.dot(m + (a-1)*mlogm,px[0])
        #eps0 = np.dot(1 + a*np.log(m[m>0]),px[0][m>0])
        
        #apprexpo = n**a
        #appreps1 = n*(1+(a-1)*np.log(n))
        #appreps0 = 1 + a*np.log(n)
        
        print("{:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(a,n,res1,appr1,appr2))
    print()



