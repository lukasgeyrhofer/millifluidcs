#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + "/..")
import growthclasses as gc


parser = argparse.ArgumentParser()
parser.add_argument("-M","--maxM",type=int,default=300)
parser.add_argument("-N","--maxN",type=float,default=10)
parser.add_argument("-n","--dN",type=float,default=.1)
parser.add_argument("-A","--maxA",type=float,default=2)
parser.add_argument("-a","--dA",type=float,default=.1)
args = parser.parse_args()


alist = np.arange(-args.maxA,args.maxA+args.dA,args.dA)
nlist = np.arange(int(np.ceil(args.maxN/args.dN))+1)*args.dN
m = np.arange(args.maxM,dtype=float)

mlogm = np.arange(args.maxM,dtype=float)
mlogm[m>0] *= np.log(m[m>0])

for n in nlist:
    px = gc.PoissonSeedingVectors(m,[n])
    for a in alist:
        ma = np.power(m[m>0],a)

        ma[ma is np.nan] = 0
        
        expo = np.dot(ma,px[0][m>0])
        eps1 = np.dot(m + (a-1)*mlogm,px[0])
        eps0 = np.dot(1 + a*np.log(m[m>0]),px[0][m>0])
        
        apprexpo = n**a
        appreps1 = n*(1+(a-1)*np.log(n))
        appreps0 = 1 + a*np.log(n)
        
        print("{:5.2f} {:5.2f} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e}".format(a,n,expo,eps1,eps0,apprexpo,appreps1,appreps0))
    print()



