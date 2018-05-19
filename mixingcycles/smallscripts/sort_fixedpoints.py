#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

import pandas as pd


def stability(ev1r,ev2r,ev1i,ev2i):
    unstable  = 0
    unstable += 1 if (ev1r**2 + ev1i**2 < 1) else 0
    unstable += 1 if (ev2r**2 + ev2i**2 < 1) else 0
    unstable += 0 if (ev1i == ev2i == 0) else 3

    return unstable

def alloweddilution(dilution,dilutionstep,threshold = 1e-10):
    if not dilutionstep is None:
        if (dilution / dilutionstep - dilution // dilutionstep)**2 < threshold:
            return True
        else:
            return False
    else:
        return True


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-o","--outbasename",default=None)
parser.add_argument("-C","--coex",default=False,action="store_true")
parser.add_argument("-D","--dilutionminstep",default=None,type=float)
args = parser.parse_args()

data = None
if len(args.infiles) == 0:  raise IOError("no input files")

for fn in args.infiles:
    try:
        newdata = np.genfromtxt(fn)
    except:
        continue

    if not data is None:
        data = np.concatenate([data,newdata])
    else:
        data = newdata[:,:]

if data is None:    raise IOError("no input data")


if args.outbasename is None:
    newname = args.infiles[0][:-2]
else:
    newname = args.outbasename
    
if not args.coex:
    data_stable          = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 2     and alloweddilution(x[0],args.dilutionminstep)])
    data_stablecomplex   = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 5     and alloweddilution(x[0],args.dilutionminstep)])
    data_unstable        = np.array([x for x in data if 0 <= stability(x[6],x[7],x[8],x[9]) < 2 and alloweddilution(x[0],args.dilutionminstep)])
    data_unstablecomplex = np.array([x for x in data if 3 <= stability(x[6],x[7],x[8],x[9]) < 5 and alloweddilution(x[0],args.dilutionminstep)])

    np.savetxt(newname + "_stable",          data_stable)
    np.savetxt(newname + "_stablecomplex",   data_stablecomplex)
    np.savetxt(newname + "_unstable",        data_unstable)
    np.savetxt(newname + "_unstablecomplex", data_unstablecomplex)
else:
    data_Cstable          = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 2     and     (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_Cstablecomplex   = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 5     and     (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_Cunstable        = np.array([x for x in data if 0 <= stability(x[6],x[7],x[8],x[9]) < 2 and     (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_Cunstablecomplex = np.array([x for x in data if 3 <= stability(x[6],x[7],x[8],x[9]) < 5 and     (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_stable           = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 2     and not (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_stablecomplex    = np.array([x for x in data if stability(x[6],x[7],x[8],x[9]) == 5     and not (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_unstable         = np.array([x for x in data if 0 <= stability(x[6],x[7],x[8],x[9]) < 2 and not (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])
    data_unstablecomplex  = np.array([x for x in data if 3 <= stability(x[6],x[7],x[8],x[9]) < 5 and not (x[3] > 0 and x[4] > 0) and alloweddilution(x[0],args.dilutionminstep)])

    np.savetxt(newname + "_stable",           data_stable)
    np.savetxt(newname + "_stablecomplex",    data_stablecomplex)
    np.savetxt(newname + "_unstable",         data_unstable)
    np.savetxt(newname + "_unstablecomplex",  data_unstablecomplex)

    np.savetxt(newname + "_Cstable",          data_Cstable)
    np.savetxt(newname + "_Cstablecomplex",   data_Cstablecomplex)
    np.savetxt(newname + "_Cunstable",        data_Cunstable)
    np.savetxt(newname + "_Cunstablecomplex", data_Cunstablecomplex)

