#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-T","--templatefile",default=None)
parser.add_argument("-m","--lowerthreshold",default=.02,type=float)
parser.add_argument("-M","--upperthreshold",default=2,type=float)
args = parser.parse_args()

try:
    fpTMP = open(args.templatefile,"r")
except:
    raise IOError,"could not open templatefile"

first = True
template = dict()
for line in fpTMP.readlines():
    if line[0] != "#":
        values = line.strip().split(',')
        if first:
            names = values
            first = False
        else:
            for i in range(len(values)):
                if not template.has_key(names[i]):
                    template[names[i]] = np.array([values[i]])
                else:
                    template[names[i]] = np.concatenate([template[names[i]],[values[i]]])

description = np.concatenate([[template['description'][i]] * int(template['droplet_number'][i]) for i in range(len(template['description']))])
growthrates = dict()

assert len(args.infiles) == len(description)

for i in range(len(args.infiles)):
    try:
        data = np.genfromtxt(args.infiles[i],delimiter = ',', names = True)
        t = data['time'] * 1e-3 # for whatever reason, the unit seems to be milli-hours
        b = data['Channel1_mean']

        if not growthrates.has_key(description[i]):
            growthrates[description[i]] = list()
        
        tEVEN = t[0::2]
        tODD  = t[1::2]
        bEVEN = b[0::2]
        bODD  = b[1::2]
        
        tEVEN = tEVEN[args.lowerthreshold < bEVEN]
        bEVEN = bEVEN[args.lowerthreshold < bEVEN]
        tODD  = tODD [args.lowerthreshold < bODD ]
        bODD  = bODD [args.lowerthreshold < bODD ]
        
        tEVEN = tEVEN[args.upperthreshold > bEVEN]
        bEVEN = bEVEN[args.upperthreshold > bEVEN]
        tODD  = tODD [args.upperthreshold > bODD ]
        bODD  = bODD [args.upperthreshold > bODD ]
        
        if len(tODD) >= 2:
            sx  = np.sum(tODD)
            sxx = np.sum(tODD*tODD)
            sy  = np.sum(np.log(bODD))
            sxy = np.sum(tODD * np.log(bODD))
            n   = len(tODD)
            grODD = (n * sxy - sx * sy)/(n*sxx - sx*sx)
        
            if np.isnan(grODD):
                print bODD,description[i]
            else:
                growthrates[description[i]].append(grODD)

        if len(tEVEN) >= 2:
            sx  = np.sum(tEVEN)
            sxx = np.sum(tEVEN*tEVEN)
            sy  = np.sum(np.log(bEVEN))
            sxy = np.sum(tEVEN * np.log(bEVEN))
            n   = len(tEVEN)
            grEVEN = (n * sxy - sx * sy)/(n*sxx - sx*sx)

            if np.isnan(grEVEN):
                print bEVEN,description[i]
            else:
                growthrates[description[i]].append(grEVEN)
        
    except:
        continue


for key in growthrates.iterkeys():
    gr = np.array(growthrates[key])
    #print gr
    print "{:15s} {:.4f} (± {:.4f}) 1/h".format(key,np.mean(gr),np.sqrt(np.std(gr)))
    h,b = np.histogram(gr,bins=10,range = (.2,.7),density = True)
    b = b[:-1] + np.diff(b)/2.
    np.savetxt(key + ".growthrates",np.transpose([b,h]))
