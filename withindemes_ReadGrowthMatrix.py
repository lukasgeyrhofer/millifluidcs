#!/usr/bin/env python3
#-*- coding: utf-8 -*-


'''
==================================
=  withindeme_ReadGrowthMatrix.py
==================================

    Reads a GrowthMatrix file generated before and
    makes its contents human-readable
    
    Required command-line parameter:
        -i GROWTHMATRIXFILE
    
    
    If not '-q'/'--quiet' is specified,
    then write all parameters that have been used
    to generate this growthmatrix to stderr.

    Output all values in growthmatrix file to
    either stdout or to 'outfile',
    if '-o OUTFILE'/'--outfile OUTFILE' is given
    
    Lukas Geyrhofer, l.geyrhofer@technion.ac.il

'''

import numpy as np
import argparse
import pickle
import sys
import growthclasses as gc

def array_to_str(x):
    return ' '.join(['{:e}'.format(float(y)) for y in x])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile",  default = None, required = True)
    parser.add_argument("-o", "--outfile", default = None)
    parser.add_argument("-q", "--quiet",   default = False, action = "store_true")
    args = parser.parse_args()

    try:
        g = pickle.load(open(args.infile,'rb'),encoding = 'bytes')
    except:
        raise IOError("Could not open pickle file")

    
    
    if not args.quiet:
        # output of all parameters
        
        fp_params = sys.stderr
        
        fp_params.write("==================================================================\n")
        fp_params.write("  growthmatrix file : {}\n".format(args.infile)                      )
        fp_params.write("  dynamics type :     {}\n".format(type(g).__name__)                 )
        fp_params.write("  GMshape :           {}\n".format(np.shape(g.growthmatrix))         )
        fp_params.write("  GMgridX :           {}\n".format(g.growthmatrixgrid[0])            )
        fp_params.write("  GMgridY :           {}\n".format(g.growthmatrixgrid[1])            )
        fp_params.write("==================================================================\n")
        fp_params.write("\n"                                                                  )

        fp_params.write(str(g))
        fp_params.write("\n")
        fp_params.write("==================================================================\n")
        fp_params.write("\n")
        fp_params.write("  command-line parameters to generate this dynamics:\n"              )


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
        
        generateString += "\n"
        
        fp_params.write(generateString)
        fp_params.write("\n")
        fp_params.write("==================================================================\n")
        
    
    # output of all values
    if not args.outfile is None:
        try:
            fp_values = open(args.outfile,"w")
        except:
            raise IOError("could not open file '{}' to write.".format(args.outfile))
    else:
        fp_values = sys.stdout
    
    if isinstance(g.growthmatrixgrid,int):
        for x in range(g.growthmatrixgrid):
            for y in range(g.growthmatrixgrid):
                fp_values.write('{} {} {} {}\n'.format(x,y,g.growthmatrix[x,y,0],g.growthmatrix[x,y,1]))
            fp_values.write('\n')
    elif isinstance(g.growthmatrixgrid,(tuple,list,np.ndarray)):
        for i,x in enumerate(g.growthmatrixgrid[0]):
            for j,y in enumerate(g.growthmatrixgrid[1]):
                fp_values.write('{} {} {} {}\n'.format(x,y,g.growthmatrix[i,j,0],g.growthmatrix[i,j,1]))
            fp_values.write('\n')
    else:
        raise NotImplementedError
    
    if not args.outfile is None:
        fp_values.close()


if __name__ == "__main__":
    main()
    
