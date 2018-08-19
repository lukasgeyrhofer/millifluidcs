#!/usr/bin/env python3

import numpy as np
import argparse
import pickle
import growthclasses as gc

def array_to_str(x):
    return ' '.join(['{:e}'.format(y) for y in x])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile",default=None)
    args = parser.parse_args()

    try:
        g = pickle.load(open(args.infile,'rb'),encoding = 'bytes')
    except:
        raise IOError("Could not open and load from pickle file")

    print("==================================================================")
    print("  growthmatrix file : {}".format(args.infile)                      )
    print("  dynamics type :     {}".format(type(g).__name__)                 )
    print("  GMshape :           {}".format(np.shape(g.growthmatrix))         )
    print("  GMgridX :           {}".format(g.growthmatrixgrid[0])            )
    print("  GMgridY :           {}".format(g.growthmatrixgrid[1])            )
    print("==================================================================")
    print(""                                                                  )

    print(g.ParameterString())
    print("")
    print("==================================================================")
    print("")
    print("  command-line parameters to generate this dynamics:"              )


    generateString = "   -G '" + type(g).__name__[14:] + "'"
    # default parameters for strains
    excludeParameters  = ['growthrates','yieldfactors','substrateconcentration','mixingtime']
    # I/O parameters from computing growthmatrix
    excludeParameters += ['infile','outfile','verbose']
    # other parameters from computing GrowthMatrix
    excludeParameters += ['GrowthDynamics','maxsize','step']
    
    if 'growthrates' in g._GrowthDynamics__kwargs_for_pickle.keys():
        generateString += " -a " + array_to_str(g.growthrates)
    if 'yieldfactors' in g._GrowthDynamics__kwargs_for_pickle.keys():
        generateString += " -y " + array_to_str(g.yieldfactors)
    if 'substrateconcentration' in g._GrowthDynamics__kwargs_for_pickle.keys():
        generateString += " -S " + str(g.env.substrate)
    if 'mixingtime' in g._GrowthDynamics__kwargs_for_pickle.keys():
        generateString += " -T " + str(g.env.mixingtime)
    
    additionalParameters = ""
    for key,values in g._GrowthDynamics__kwargs_for_pickle.items():
        if not key in excludeParameters:
            if isinstance(values,(list,np.ndarray)):
                additionalParameters += " {} ".format(key) + array_to_str(values)
            else:
                additionalParameters += " {} {}".format(key,values)
    if len(additionalParameters) > 0:
        generateString += " -P " + additionalParameters
    
    print(generateString)
    print("")
    print("==================================================================")


if __name__ == "__main__":
    main()
    
