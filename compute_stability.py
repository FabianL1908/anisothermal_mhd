#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:58:52 2020

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:58:29 2019

@author: boulle
"""

import sys
import os
import gc

from firedrake import *
import numpy as np
from defcon import *
from defcon import backend, StabilityTask
from defcon.cli.common import fetch_bifurcation_problem
from petsc4py import PETSc

RB = __import__("linear_eigenvalue")

# Branch
comm = COMM_WORLD

# Construct mono-3d problem
problem = RB.EVRayleighBenardProblem()
mesh = problem.mesh(comm=comm)
Z = problem.function_space(mesh)

# Set-up io function
io = problem.io("output")
params = problem.parameters()
io.setup(params, problem.functionals(), Z)

# set the petsc options from the solver_parameters
solver_parameters = problem.solver_parameters(problem.parameters(), StabilityTask)
# PETSC options
opts = PETSc.Options()
for k in solver_parameters:
    opts[k] = solver_parameters[k]

Ra = 20000
Ra = 10**5
Pr = 1.0
#consts = [Ra, Pr]
consts = [max(y[0], y[-1]) for x, y in problem.parameter_values().items()]

#consts = problem.parameters()
print(consts)
#solution = io.fetch_solutions(consts, [0])[0]
solution = Function(Z)
d = problem.compute_stability(consts, 0, solution)
io.save_stability(d["stable"], d.get("eigenvalues", []), d.get("eigenfunctions", []), consts, 0)
# pvd = File("eigenfunctions/eigenfunctions.pvd")
# for e in d.get("eigenfunctions", []):
#     problem.save_pvd(e, pvd)
eigs = d.get("eigenvalues", [])
evals = list(map(complex, eigs))
evalsR = np.array([l.real for l in evals])
eigsPos = evalsR[evalsR>=0]
print(eigsPos)
print(len(eigsPos))
