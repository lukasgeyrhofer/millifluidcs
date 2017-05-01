#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import os

class DropletData(object):
    def __init__(self, templatefile = None, infiles = None, splitBackForthTrajectories = False, datacolumns = ['time','Channel1_mean'], snakelikeloading = True):
        
        if not templatefile is None:
            try:
                fptemp = open(templatefile,"r")
            except:
                raise IOError,"could not open file"
            first = True
            for line in fptemp.readlines():
                if line[0] != "#":
                    values = line.strip().split(',')
                    if first:
                        names = values
                        try:
                            IDdescription    = names.index("description")
                            IDdroplet_number = names.index("droplet_number")
                        except:
                            raise ValueError, "templatefile does not contain columns 'description' and 'droplet_number'"

                        self.__droplettype = None
                        typesinrow         = None
                        lastrow            = ""
                        first              = False
                    else:
                        if line[0] != lastrow:
                            if snakelikeloading:
                                # check if forward
                                direction      = 2 * (ord(line[0])%2) - 1   # ord('A') == 65
                                                                            # seems quite a hack, not sure how general this is
                            else:
                                direction      = 1
                            self.__droplettype = self.concat(self.__droplettype,typesinrow[::direction])
                            typesinrow         = None
                        
                        typesinrow = self.concat(typesinrow,np.repeat(values[IDdescription],int(values[IDdroplet_number])))
                        lastrow = line[0]
        
            fptemp.close()
        else:
            self.__droplettype = np.repeat("default",len(infiles),dtype=str)
        
        self.__data        = dict()
        self.__listoftypes = list(set(self.__droplettype))
        self.__datacolumns = datacolumns
            
        for filename in infiles:
            try:
                tmpdata = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
            except:
                raise IOError
            dropletID = self.filename_to_dropletID(filename)
            self.add_trajectory(dropletID, tmpdata, self.__datacolumns, splitBackForthTrajectories)
    
    
    def concat(list1,list2):
        if (list1 is None):
            return list2
        elif (list2 is None):
            return list1
        elif (list1 is None) and (list2 is None):
            return None
        else:
            return np.concatenate([list1,list2])

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
                newdatablock0 = self.concat(newdatablock0,np.array([data[column][0::2]]))
                newdatablock1 = self.concat(newdatablock1,np.array([data[column][1::2]]))                                                    
            newdatablock0 = np.transpose(newdatablock0)
            newdatablock1 = np.transpose(newdatablock1)
            self.__data[dropletLabel].append(newdatablock0)
            self.__data[dropletLabel].append(newdatablock1)
        else:
            newdatablock0 = None
            for column in columns:
                newdatablock0 = self.concat(newdatablock0,np.array([data[column][0::2]]))
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
    
    
    def CountTrajectories(self):
        r = dict()
        for key in self.__listoftypes:
            r[key] = len(self.__data[key])
        return r
    
    
    def __getattr__(self,key):
        if key == "datacolumns":
            return self.__datacolumns
        else:
            super(DropletData,self).__getattr__(key)

# working example of loading files and printing them again
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*")
    parser.add_argument("-T","--templatefile",default=None)
    parser.add_argument("-S","--splitBackForthTrajectories",default=False,action="store_true")
    args = parser.parse_args()
    
    # define data object
    data = DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = args.splitBackForthTrajectories)

    # in general, data is obtained by iterating over sets of experimental labels with their trajectories
    for dropletLabel, Trajectories in data:
        print dropletLabel, len(Trajectories)
        for trajectory in Trajectories:
                print trajectory
                
    # should also work like that:
    # given that the label 'SBW25mC-10' is defined in the templatefile
    # this is the same as the 'Trajectories' in the loop above
    print data['SBW25mC-10']


# run through working example, if program is executed directly
# don't run, if just loaded for class definitions
if __name__ == "__main__":
    main()






