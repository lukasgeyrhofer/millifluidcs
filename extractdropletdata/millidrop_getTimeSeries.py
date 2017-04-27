#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-L","--forcelen",type=int,default=None)
parser.add_argument("-S","--splitEvenOdd",default=False,action="store_true")
parser.add_argument("-T","--templatefile",default=None)
args = parser.parse_args()


if not args.templatefile is None:
    try:
        fpTMP = open(args.templatefile,"r")

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

    except:
        description = None


i = 0
for filename in args.infiles:
    try:
        data = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
        if args.splitEvenOdd:
            newdata0 = np.transpose(np.array([data['time'][0::2],data['Channel1_mean'][0::2]]))
            newdata1 = np.transpose(np.array([data['time'][1::2],data['Channel1_mean'][1::2]]))
            
            if not description is None:
                newfilename0 = filename.split(".")[0] + '-' + description[i] + '-0.data'
                newfilename1 = filename.split(".")[0] + '-' + description[i] + '-1.data'
                i += 1
            else:
                newfilename0 = filename.replace(".csv","-0.data")
                newfilename1 = filename.replace(".csv","-1.data")
            
            np.savetxt(newfilename0,newdata0,fmt='%.6e')
            np.savetxt(newfilename1,newdata1,fmt='%.6e')
        else:
            newdata = np.transpose(np.array([data['time'],data['Channel1_mean']]))
            if not args.forcelen is None:
                if len(newdata) < args.forcelen:
                    zeros = np.zeros(shape = (args.forcelen- len(newdata),2))
                    newdata = np.concatenate((newdata,zeros))
                    print "forcelength:",filename
                    
            if not description is None:
                newfilename = filename.split(".")[0] + '-' + description[i] + '.data'
                i += 1
            else:
                newfilename = filename.replace("csv","data")
            np.savetxt(newfilename,newdata,fmt='%.6e')
    except:
        continue
