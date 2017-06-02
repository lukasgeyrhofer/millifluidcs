#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import os

class DropletData(object):
    def __init__(self, templatefile = None, infiles = None, splitBackForthTrajectories = False, datacolumns = ['time','Channel1_mean'], snakelikeloading = True, restrictionfile = None):
        
        if infiles is None:
            raise IOError("datafiles required")
        
        if not templatefile is None:
            self.load_templatefile(templatefile,snakelikeloading)
        else:
            # if no templatefile, group everything under label 'default' with 'unkown' well ID
            self.__droplettype = np.repeat("default",len(infiles))
            self.__well        = np.repeat("unknown",len(infiles))
        
        # ===============================================================
        # = variable initialization
        # ===============================================================
        self.__listoftypes         = list(set(self.__droplettype))  # list of all possible "experiments", ie. labels of different wells
        self.__data                = dict()                         # store all trajectories, each keys of the dictionary are the experiment labels, "values" of such an entry is list of np-arrays
        self.__datacolumns         = datacolumns                    # keep only those columns from the original datafile
        self.__datarestrictions    = list()                         # list of all restrictions
        self.__permittedoperations = {"min":2,"max":2,"end":2,"start":2,"exclude":0}
                                                                    # hardcoded allowed restriction types: min/max chop trajectories if the fall below or rise above the threshold value, end/start is for time


        # ===============================================================
        # = iterate over all separate droplet files
        # ===============================================================
        for filename in infiles:
            dropletID = self.filename_to_dropletID(filename)
            try:
                tmpdata = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
            except:
                print("Error while loading file '{:s}'. Continuing ...".format(filename))
                continue
            self.add_trajectory(dropletID, tmpdata, self.__datacolumns, splitBackForthTrajectories)
        
        
        # ===============================================================
        # = restrictions on data
        # ===============================================================
        
        if not restrictionfile is None:
            self.load_restrictions_from_file(restrictionfile)
        



    # ===============================================================
    # = generate list of all types of experiments from templatefile =
    # ===============================================================
    def load_templatefile(self,templatefile = None, snakelikeloading = True):
        try:
            fptemp = open(templatefile,"r")
        except:
            raise IOError("could not open file")
        first = True
        for line in fptemp.readlines():
            if line[0] != "#":
                values = line.strip().split(',')
                if first:
                    names = values
                    try:
                        IDdescription    = names.index("description")
                        IDdroplet_number = names.index("droplet_number")
                        IDwell           = names.index("well")
                    except:
                        raise ValueError("templatefile does not contain columns 'description' and 'droplet_number'")

                    self.__droplettype = None
                    self.__well        = None
                    typesinrow         = None
                    wellsinrow         = None
                    lastrow            = '.'
                    first              = False
                else:
                    if line[0] != lastrow:
                        if snakelikeloading:
                            # check if forward
                            direction      = 2 * (ord(lastrow)%2) - 1   # ord('A') == 65
                                                                        # seems quite a hack, not sure how general this is
                        else:
                            direction      = 1

                        self.__droplettype = self.concat(self.__droplettype,typesinrow,direction2 = direction)
                        typesinrow         = None
                        self.__well        = self.concat(self.__well,wellsinrow,direction2 = direction)
                        wellsinrow         = None
                    
                    typesinrow = self.concat(typesinrow,np.repeat(values[IDdescription],int(values[IDdroplet_number])))
                    wellsinrow = self.concat(wellsinrow,np.repeat(values[IDwell],int(values[IDdroplet_number])))
                    lastrow = line[0]
        
        # add the last buffer 'typesinrow'
        if snakelikeloading:    direction      = 2 * (ord(lastrow)%2) - 1
        else:                   direction      = 1
        self.__droplettype = self.concat(self.__droplettype,typesinrow,direction2 = direction)
        self.__well        = self.concat(self.__well,wellsinrow,direction2 = direction)
        
        # close templatefile
        fptemp.close()

    # ===============================================================
    # = load data from a single dropletfile and add to internal datastructure
    # ===============================================================
    def add_trajectory(self,dropletID = None, data = None, columns = None, splitBackForthTrajectories = False):
        if data is None:
            raise ValueError("need timestamps for experimental data")

        for column in columns:
            if not column in data.dtype.names:
                raise ValueError("No field of name '{}'. Possible values are ('".format(column) + "', '".join(data.dtype.names) + "')")
        dropletLabel = self.dropletID_to_label(dropletID)
        if not self.__data.has_key(dropletLabel):
            self.__data[dropletLabel] = list()
        
        if splitBackForthTrajectories:        
            newdatablock0 = dict()
            newdatablock1 = dict()
            for column in columns:
                newdatablock0[column] = np.array(data[column][0::2])
                newdatablock1[column] = np.array(data[column][1::2])                                                    
            self.__data[dropletLabel].append(newdatablock0)
            self.__data[dropletLabel].append(newdatablock1)
        else:
            newdatablock0 = dict()
            for column in columns:
                newdatablock0[column] = np.array(data[column])
            self.__data[dropletLabel].append(newdatablock0)
            
    # ===============================================================
    # = helper routines
    # ===============================================================

    def concat(self,list1 = None,list2 = None, direction1 = 1, direction2 = 1):
        if (list1 is None) and (list2 is None):
            return None
        elif (list1 is None):
            return list2[::direction2]
        elif (list2 is None):
            return list1[::direction1]
        else:
            return np.concatenate([list1[::direction1],list2[::direction2]])


    def filename_to_dropletID(self,filename):
        # in current implementation, filename is just 'dropletnumber.CSV'
        # change here if necessary
        r = int(os.path.basename(filename).split('.')[0])
        return r


    def dropletID_to_label(self,dropletID = None):
        if dropletID is None:
            raise KeyError("dropletID is None")
        try:
            return self.__droplettype[dropletID]
        except:
            raise KeyError("did not find experiment type for droplet '{}'".format(dropletID))
    
    
    # ===============================================================
    # = output of specific information from internal datastructure
    # ===============================================================

    def CountTrajectories(self):
        r = dict()
        for key in self.__listoftypes:
            r[key] = len(self.__data[key])
        return r
    
    
    def getDropletTypes(self):
        return self.__droplettype
    
    
    def getChannelIndex(self,key):
        if key in self.__datacolumns:
            return self.__datacolumns.index(key)
        else:
            return None
    
    
    # ===============================================================
    # routines to work with restricted data (ie. a lower cutoff for the droplet signal, or a maximal time)
    # ===============================================================
    
    def set_restriction(self,restictiontype, applies_to = None, params = None):
        if restrictiontype in self.__permittedoperations.keys():
            if applies_to is None:
                applies_to = "all"
            if (applies_to in self.__droplettype) or (applies_to in set(self.__well)) or (applies_to == "all"):
                if len(params) == self.__permittedoperations[restrictiontype]:
                    newrestriction = list(restrictiontype,applies_to,params)
                    self.__datarestrictions.append(newrestriction)
        else:
            raise ValueError("cannot apply restriction to data (column: '{:s}', operation: {:s})".format(column,operation))


    def remove_all_restrictions(self):
        self.__datarestrictions = list()
        
        
    def load_restrictions_from_file(self,filename):
        try:
            fp = open(filename,"r")
        except:
            raise IOError("could not open file '{:s}' to load restrictions".format(filename))
        for line in fp.readlines():
            values = line.split()
            if len(values) >= 2:
                if values[0] in self.__permittedoperations.keys():
                    if (values[1] in self.__droplettype) or (values[1] in set(self.__well)) or (values[1] == "all"):
                        if len(values) == 2 + self.__permittedoperations[values[0]]:
                            self.set_restriction(values[0],values[1],values[2:])
        fp.close()
    
    
    def write_restrictions_to_file(self,filename = None):
        try:
            if filename is None:
                fp = sys.stdout
            else:
                fp = open(filename,"w")
        except:
            raise IOError("could not open file '{:s}' to save restrictions".format(filename))
        for restriction in self.__datarestrictions:
            print >> fp, restriction[0] + " " + restriction[0] + " " + " ".join([str(x) for x in restriction[2]])
        if not filename is None:
            fp.close()

    # ===============================================================
    # all output of data is funneled through this routine
    # only return datapoints that match all criteria in the restrictionfile
    # ===============================================================
    
    def pattern_and(self,pattern1,pattern2):
        if pattern1 is None:
            return pattern2
        else:
            return np.logical_and(pattern1,pattern2)
    
    def restricted_data(self,key):
        r = list()
        for datablock in self.__data[key]:
            rdata = None
            for column in self.__datacolumns:
                if rdata is None:
                    rdata = np.array([datablock[column]])
                else:
                    rdata = np.concatenate([rdata,np.array([datablock[column]])],axis = 0)
            rdata = np.transpose(rdata)
            pattern = None
            if len(self.__datarestrictions) > 0:
                for restriction in self.__datarestrictions:
                    if   restriction[0] == "max":
                        pattern = self.pattern_and(pattern, datablock[restriction[2,0]] < float(restriction[2,1]) )
                    elif restriction[0] == "min":
                        pattern = self.pattern_and(pattern, datablock[restriction[2,0]] > float(restriction[2,1]) )
                    elif restriction[0] == "end":
                        pattern = self.pattern_and(pattern, datablock[restriction[2,0]] < (1-restriction[2,1]) * datablock[restriction[2,0]][-1] )
                    elif restriction[0] == "start":
                        pattern = self.pattern_and(pattern, datablock[restriction[2,0]] > (1+restriction[2,1]) * datablock[restriction[2,0]][0] )
                    elif restriction[0] == "exclude":
                        continue
                        #pattern = self.pattern_and(pattern, )
                
                #for restriction in self.__datarestrictions:
                    #if len(restriction) == 4:
                        #if restriction[3] == key:
                            #continue
                    #IDrestiction = self.__datacolumns.index(restriction[0])
                    #if restriction[2] == "max":
                        #if pattern is None:
                            #pattern = (datablock[restriction[0]] < restriction[1])
                        #else:
                            #pattern = np.logical_and(pattern, datablock[restriction[0]] < restriction[1])
                    #elif restriction[2] == "min":
                        #if pattern is None:
                            #pattern = (datablock[restriction[0]] > restriction[1])
                        #else:
                            #pattern = np.logical_and(pattern,datablock[restriction[0]] > restriction[1])
                    #elif restriction[2] == "end":
                        #if pattern is None:
                            #pattern = (datablock[restriction[0]] < (1-restriction[1])* datablock[restriction[0]][-1])
                        #else:
                            #pattern = np.logical_and(pattern,datablock[restriction[0]] < (1-restriction[1])* datablock[restriction[0]][-1])
                    #elif restriction[2] == "start":
                        #if pattern is None:
                            #pattern = (datablock[restriction[0]] > (1+restriction[1]) * datablock[restriction[0]][0])
                        #else:
                            #pattern = np.logical_and(pattern,datablock[restriction[0]] (1+restriction[1]) * datablock[restriction[0]][0])
                        
                pattern = np.transpose(np.repeat([pattern],len(rdata[0]),axis = 0))
                rdata   = np.reshape(rdata[pattern],(len(rdata[pattern])/len(self.__datacolumns),len(self.__datacolumns)))
            r.append(rdata)
        return r
          
    
    # ===============================================================
    # = object logic                                                =
    # ===============================================================
    
    # after object initialization, this is how the dataobject reacts
    # example code lines start with '>'
    # > data = DropletData()
    
    # get all trajectories of a certain experiment as list,
    # when EXPERIMENTLABEL is a string also contained in TEMPLATEFILE as label for a well
    # > data[EXPERIMENTLABEL]
    def __getitem__(self,key):
        if key in self.__listoftypes:
            return self.restricted_data(key)
    
    # iterate over the whole dataset
    # > for ExperimentLabel, Trajectories in data:
    # >     print ExperimentLabel
    # >     for trajectory in Trajectories:
    # >         print trajectory
    def __iter__(self):
        for key in self.__listoftypes:
            yield key,self[key]

    # catch any specific attributes of the DropletData object
    # so far only single attribute implemented
    # > print data.datacolumns
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






