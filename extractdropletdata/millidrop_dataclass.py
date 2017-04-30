#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import os

class DropletData(object):
    def __init__(self, templatefile = None, infiles = None, splitBackForthTrajectories = False, datacolumns = ['time','Channel1_mean']):
        
        if not templatefile is None:
            try:
                fptemp = open(templatefile,"r")
            except:
                raise IOError,"could not open file"
        
            self.__droplettype = None
            first = True
            for line in fptemp.readlines():
                if line[0] != "#":
                    values = line.strip().split(',')
                    if first:
                        names = values
                        first = False
                        try:
                            IDdescription    = names.index("description")
                            IDdroplet_number = names.index("droplet_number")
                        except:
                            raise ValueError, "templatefile does not contain columns 'description' and 'droplet_number'"
                    else:
                        if self.__droplettype is None:
                            self.__droplettype = np.repeat(values[IDdescription],int(values[IDdroplet_number]))
                        else:
                            self.__droplettype = np.concatenate([self.__droplettype,np.repeat(values[IDdescription],int(values[IDdroplet_number]))])
        
            fptemp.close()
        else:
            self.__droplettype = np.repeat("default",len(infiles),dtype=str)
        
        self.__data        = dict()
        self.__listoftypes = list(set(self.__droplettype))
            
        for filename in infiles:
            try:
                tmpdata = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
            except:
                raise IOError
            dropletID = self.filename_to_dropletID(filename)
            self.add_trajectory(dropletID, tmpdata, datacolumns, splitBackForthTrajectories)
            

    def filename_to_dropletID(self,filename):
        # in current implementation, filename is just 'dropletnumber.CSV'
        # change here if necessary
        r = int(os.path.basename(filename).split('.')[0])
        return r

    def dropletID_to_label(self,dropletID = None):
        if dropletID is None:
            raise KeyError
        try:
            return self.__droplettype[dropletID]
        except:
            raise KeyError


    def add_trajectory(self,dropletID = None, data = None, columns = None, splitBackForthTrajectories = False):
        
        if data is None:
            raise ValueError,"need timestamps for experimental data"

        dropletLabel = self.dropletID_to_label(dropletID)
        if not self.__data.has_key(dropletLabel):
            self.__data[dropletLabel] = list()
        
        if splitBackForthTrajectories:        
            newdatablock0 = None
            newdatablock1 = None
            for column in columns:
                if newdatablock0 is None:
                    newdatablock0 = np.array([data[column][0::2]])
                    newdatablock1 = np.array([data[column][1::2]])
                else:
                    newdatablock0 = np.concatenate([newdatablock0,np.array([data[column][0::2]])])
                    newdatablock1 = np.concatenate([newdatablock1,np.array([data[column][1::2]])])                                                    
            newdatablock0 = np.transpose(newdatablock0)
            newdatablock1 = np.transpose(newdatablock1)
            self.__data[dropletLabel].append(newdatablock0)
            self.__data[dropletLabel].append(newdatablock1)
        else:
            newdatablock0 = None
            for column in columns:
                if newdatablock0 is None:
                    newdatablock0 = np.array([data[column][0::2]])
                else:
                    newdatablock0 = np.concatenate([newdatablock0,np.array([data[column][0::2]])])
            newdatablock0 = np.transpose(newdatablock0)
            self.__data[dropletLabel].append(newdatablock0)
            
            
    def __getitem__(self,key):
        if key in self.__listoftypes:
            return self.__data[key]
        else:
            super(DropletData,self).__getitem__(key)
    
    def __iter__(self):
        for key in self.__listoftypes:
            yield key,self[key]


# working example of loading files

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*")
    parser.add_argument("-T","--templatefile",default=None)
    parser.add_argument("-S","--splitBackForthTrajectories",default=False,action="store_true")
    args = parser.parse_args()
    
    data = DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = args.splitBackForthTrajectories)

    for dropletLabel, Trajectories in data:
        print dropletLabel, len(Trajectories)
        for trajectory in Trajectories:
                print trajectory
    
if __name__ == "__main__":
    main()






