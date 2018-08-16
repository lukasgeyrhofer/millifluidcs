#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import growthclasses as gc

def main():
    parser = argparse.ArgumentParser()
    parser = gc.AddGrowthParameters(parser)

    parser_gc = parser.add_argument_group("==== GrowthDynamics ====")
    parser_gc.add_argument("-d","--GrowthDynamics",default="")
    parser_gc.add_argument("-L","--ParameterList",nargs="*",default=[])

    parser_alg = parser.add_argument_group(description = "==== Parameters for algorithm ====")
    parser_alg.add_argument("-t","--TimeIntegratorStep",default=1e-3,type=float)
    parser_alg.add_argument("-O","--TimeIntegratorOutput",default=10,type=int)
    parser_alg.add_argument("-M","--IntegrationMethod",choices = ['ownRK4','SciPy'],default = 'ownRK4')

    parser_ic = parser.add_argument_group(description = "==== Initial conditions ====")
    parser_ic.add_argument("-N","--initialconditions",default=[1],type=float,nargs="*")

    args = parser.parse_args()

    g = gc.AssignGrowthDynamics(**vars(args))

    # compute trajectory
    traj = g.Trajectory(args.initialconditions,TimeOutput=True)
    # output
    for x in traj:
        print ' '.join(['{:14.6e}'.format(y) for y in x])


if __name__ == "__main__":
    main()
