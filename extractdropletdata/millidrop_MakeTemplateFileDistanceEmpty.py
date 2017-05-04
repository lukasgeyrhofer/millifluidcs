#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

#import millidrop_dataclass as mdc


def SnakeLikeLoading(template):
    ret          = list()
    currentblock = list()
    lastwell     = ""
    for currentitem in template:
        if (lastwell != currentitem[0][0]) and (len(currentblock) > 0):
            direction = 1 - 2 * (ord(currentitem[0][0])%2)
            currentblock = currentblock[::direction]
            ret += currentblock
            currentblock = list()
        currentblock.append(currentitem)
        lastwell = currentitem[0][0]
    ret += currentblock
    return ret


def is_empty(currentindex):
    if template[currentindex][labelID] == args.emptyType:   return True
    else:                                                   return False


def findNextNonEmpty(currentindex):
    if not is_empty(currentindex):
        return currentindex,0
    else:
        j                 = 0
        emptydropletcount = 0
        while True:
            if currentindex+j < len(template):
                if is_empty(currentindex + j):
                    emptydropletcount += int(template[currentindex+j][dnID])
                    j += 1
                else:
                    return currentindex + j,emptydropletcount
            else:
                return None,emptydropletcount
    


parser = argparse.ArgumentParser()
ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
ioparser.add_argument("-i","--infiles",nargs="*")
ioparser.add_argument("-t","--templatefile",default=None)
ioparser.add_argument("-o","--outfile",default=None)

parser.add_argument("-E","--emptyType",default="KBIPTG",type = str)
args = parser.parse_args()


global template
template    = list()
maxdistance = 10
names       = None

fp = open(args.templatefile)
for line in fp.readlines():
    values = line.split(',')
    if line[0] != "#":
        if names is None:
            global labelID
            global dnID
            names   = [v.strip() for v in values]
            labelID = names.index("description")
            dnID    = names.index("droplet_number")
        else:
            template.append([v.strip() for v in values])
fp.close()

template    = SnakeLikeLoading(template)
newtemplate = list([names])
lastlabel   = None

i = 0
while i < len(template):
    if not is_empty(i):
        newtemplate += [template[i]]
        lastlabel    = template[i][labelID]
        i           += 1
    else:
        ne,totaldroplets = findNextNonEmpty(i)
        print i,ne,totaldroplets
        #exit(1)
        if not ne is None:
            nextlabel = template[ne][labelID]
            for j in range(totaldroplets/2):
                curitem          = template[i][:]
                curitem[labelID] = "empty-{:02d}-{:s}".format(j+1,lastlabel)
                curitem[dnID]    = '1'
                newtemplate     += [curitem]
            for j in range(totaldroplets - totaldroplets/2):
                curitem          = template[i][:]
                curitem[labelID] = "empty-{:02d}-{:s}".format(totaldroplets/2 - j,nextlabel)
                curitem[dnID]    = '1'
                newtemplate     += [curitem]
            i = ne
        else:
            for j in range(totaldroplets):
                curitem          = template[i][:]
                curitem[labelID] = "empty-{:02d}-{:s}".format(j,lastlabel)
                curitem[dnID]    = '1'
                newtemplate     += [curitem]
            i = len(template)


newtemplate = SnakeLikeLoading(newtemplate)
if args.outfile is None:    np.savetxt(sys.stdout,   newtemplate, delimiter=',', fmt = '%s')
else:                       np.savetxt(args.outfile, newtemplate, delimiter=',', fmt = '%s')


            
