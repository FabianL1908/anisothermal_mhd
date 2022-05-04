#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:20:44 2020

@author: boulle
"""

from firedrake import *
from defcon import *
import numpy as np

RB = __import__("linear_eigenvalue")

# Branch
comm = COMM_WORLD

PlotGuesses = True

# Construct mono-3d problem
problem = RB.EVRayleighBenardProblem()
mesh = problem.mesh(comm=comm)
Z = problem.function_space(mesh)

# Set-up io function
io = problem.io("output")
params = problem.parameters()
io.setup(params, problem.functionals(), Z)

# Get eigenvectors
#consts = [Ra, Pr]
consts = [max(y[0], y[-1]) for x, y in problem.parameter_values().items()]
Ra = consts[0]
Pr = consts[1]

stabs = io.fetch_stability(consts, [0], fetch_eigenfunctions=True)
stab = stabs[0]
eigenfunctions = stab["eigenfunctions"]

init_guess = []
x = SpatialCoordinate(mesh)
u0, p0, T0, B0, E0 = problem.equilibrium_solution()
trivial = Function(Z)
trivial.sub(0).interpolate(u0)
trivial.sub(1).interpolate(p0)
trivial.sub(2).interpolate(T0)
trivial.sub(3).interpolate(B0)
trivial.sub(4).interpolate(E0)

#trivial = project(as_vector([0, 0, x[1]*(2-x[1])*Ra*Pr/2, 1-x[1], 0.0, 1.0, 0.0]), Z)
(u_t, p_t, T_t, B_t, E_t) = split(trivial)
io = problem.io("initial_guess")
params = problem.parameters()
io.setup(params, problem.functionals(), Z)

if PlotGuesses:
    pvd = File("initial_guess/solution/guess.pvd")

# Trivial state
i = 0
init = Function(Z)
io.save_solution(init, [], consts, i)
init.assign(trivial)
io.save_solution(init, [], consts, i)
if PlotGuesses:
    problem.save_pvd(init, pvd, params)

eigval = stab["eigenvalues"]
i = 1
i_eig = 0
for i in range(len(eigenfunctions)):
    eigenfunction = eigenfunctions[i]
#    eigenfunction.split()[3].interpolate(Constant((0.0, 0.0)))
#    eigenfunction.split()[4].interpolate(Constant(0.0))
    init = Function(Z)
    (u, p, T, B, E) = split(eigenfunction)
    # Rescale the state
    if norm(u) >= 1e-10 and eigval[i].real > 0:
        print(eigval[i])
#        norm_factor = min(norm(T_t)/norm(T), norm(B_t)/norm(B))
        norm_T = float(norm(T_t)/norm(T))
        norm_B = float(norm(B_t)/norm(B))
        eigenfunction.split()[0].interpolate(norm_T * eigenfunction.split()[0])
        eigenfunction.split()[1].interpolate(norm_T * eigenfunction.split()[1])
        eigenfunction.split()[2].interpolate(norm_T * eigenfunction.split()[2])
        eigenfunction.split()[3].interpolate(norm_B * eigenfunction.split()[3])
        eigenfunction.split()[4].interpolate(Constant(0.0))
        init.assign(eigenfunction+trivial)
        init_guess.append(u)
        io.save_solution(init, [], consts, i_eig)
        if PlotGuesses:
#            problem.save_pvd(eigenfunction, pvd, params)
            problem.save_pvd(init, pvd, params)
        i_eig += 1
    i += 1
print(len(init_guess))
