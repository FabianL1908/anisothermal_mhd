from __future__ import absolute_import, print_function

import sys
import os
import gc

from firedrake import *
from firedrake.mg.utils import get_level
import numpy as np
from defcon import *
from defcon import backend, StabilityTask
from defcon.cli.common import fetch_bifurcation_problem
from petsc4py import PETSc

import json
import glob
import shutil

import csv
import matplotlib.pyplot as plt

RB = __import__("rayleigh-benard")

# Branch

#branchids = [44]
#branchids = [34,54,63,64,89]
branchids = [0, 38, 2, 42, 166, 161]

comm = COMM_WORLD

## Construct mono-3d problem
problem = RB.RayleighBenardProblem()
mesh = problem.mesh(comm=comm)
Z = problem.function_space(mesh)
functionals = problem.functionals()

# Set-up io function
io = problem.io("output")
params = problem.parameters()
io.setup(params, problem.functionals(), Z)

try:
    os.makedirs("diagram_u")
    os.makedirs("diagram_T")
    os.makedirs("diagram_B")
    print("Directory Created ")
except FileExistsError:
    print("Directory already exists")

for branchid in branchids:
    # Parameters
    knownparams = io.known_parameters(fixed={"Pr":1, "S": 100, "Pm": 10}, branchid=branchid)
    knownparams_Ra = np.array([l[0] for l in knownparams])
    Nu = np.array([])
    NT = np.array([])
    NB = np.array([])
    for param in knownparams:
        print(param)
        print("Computing functional for parameters %s, branchid = %d" % (str(param[0]), branchid), flush=True)
        solution = io.fetch_solutions(param, [branchid])[0]
        
        funcs = []
        for functional in functionals:
            func = functional[0]
            j = func(solution, param)
            funcs.append(j)
        
        Nu = np.append(Nu, funcs[0])
        NT = np.append(NT, funcs[1])
        NB = np.append(NB, funcs[2])
    
    # save to text file
    knownparams_Ra = knownparams_Ra.reshape((len(knownparams_Ra), 1))
    Nu = Nu.reshape((len(Nu), 1))
    NT = NT.reshape((len(NT), 1))
    NB = NB.reshape((len(NB), 1))
    np.savetxt("diagram_u/%d.csv"%branchid, np.hstack((knownparams_Ra, Nu)), delimiter=",")
    np.savetxt("diagram_T/%d.csv"%branchid, np.hstack((knownparams_Ra, NT)), delimiter=",")
    np.savetxt("diagram_B/%d.csv"%branchid, np.hstack((knownparams_Ra, NB)), delimiter=",")



if True:
    for func_idx, dgrm_type in enumerate(["u", "T", "B"]):
        for branchid in branchids:
            with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r') as f:
                data = list(csv.reader(f, delimiter=","))
            data = np.array(data)
            data = data.astype(np.float32)
            data = data.T
            plt.plot(data[0], data[1])
        plt.xlabel(r"$\mathrm{Ra}$")
        plt.ylabel(functionals[func_idx][2])
        plt.savefig(f'diagram_{dgrm_type}.png', dpi=400)
    
#shutil.make_archive("/home/boulle/Documents/diagram_data", 'zip', "diagram_data")


