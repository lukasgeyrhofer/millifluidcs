#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import pandas as pd


def stability(ev1r,ev2r,ev1i,ev2i):
    unstable  = 0
    unstable += (ev1r**2 + ev1i**2 < 1)
    unstable += (ev2r**2 + ev2i**2 < 1)
    unstable += 0 if (ev1i == ev2i == 0) else 3

    return unstable

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-C","--coex",default=False,action="store_true")
args = parser.parse_args()

data = None

for fn in args.infiles:
    try:
        newdata = np.gen_from_txt(fn)
    except:
        continue

    if not data is None:
        data = np.concatenate([data,newdata])
    else:
        data = newdata[:,:]



if not args.coex:
    data_stable          = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 2])
    data_stablecomplex   = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 5])
    data_unstable        = np.array([x for x in data if 0 <= stability(x[6],x[7],x[8],x[9])**2 < 2])
    data_unstablecomplex = np.array([x for x in data if 3 <= stability(x[6],x[7],x[8],x[9])**2 < 5])

    newname = args.infiles[:-2]

    np.save_txt(data_stable,newname + "_stable")
    np.save_txt(data_stablecomplex,newname + "_stablecomplex")
    np.save_txt(data_unstable,newname + "_unstable")
    np.save_txt(data_unstablecomplex,newname + "_unstablecomplex")
else:
    data_Cstable         = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 2 and (x[3] > 0 and x[4] > 0)])
    data_Cstablecomplex   = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 5 and (x[3] > 0 and x[4] > 0)])
    data_Cunstable        = np.array([x for x in data if 0 <= stability(x[6],x[7],x[8],x[9])**2 < 2 and (x[3] > 0 and x[4] > 0)])
    data_Cunstablecomplex = np.array([x for x in data if 3 <= stability(x[6],x[7],x[8],x[9])**2 < 5 and (x[3] > 0 and x[4] > 0)])
    data_stable          = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 2 and not (x[3] > 0 and x[4] > 0)])
    data_stablecomplex   = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 5 and not (x[3] > 0 and x[4] > 0)])
    data_unstable        = np.array([x for x in data if 0 <= stability(x[6],x[7],x[8],x[9])**2 < 2 and not (x[3] > 0 and x[4] > 0)])
    data_unstablecomplex = np.array([x for x in data if 3 <= stability(x[6],x[7],x[8],x[9])**2 < 5 and not (x[3] > 0 and x[4] > 0)])
    
    newname = args.infiles[:-2]

    np.save_txt(data_stable,newname + "_stable")
    np.save_txt(data_stablecomplex,newname + "_stablecomplex")
    np.save_txt(data_unstable,newname + "_unstable")
    np.save_txt(data_unstablecomplex,newname + "_unstablecomplex")

    np.save_txt(data_Cstable,newname + "_Cstable")
    np.save_txt(data_Cstablecomplex,newname + "_Cstablecomplex")
    np.save_txt(data_Cunstable,newname + "_Cunstable")
    np.save_txt(data_Cunstablecomplex,newname + "_Cunstablecomplex")


