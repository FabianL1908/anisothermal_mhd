#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import gc
from multiprocessing import Pool

from firedrake import *
import numpy as np
from defcon import *
from defcon import backend, StabilityTask
from defcon.cli.common import fetch_bifurcation_problem
from petsc4py import PETSc
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

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

Ras = linspace(10**3, 10**5, 400)
Pr = 1.0
S = 1.0
Pm = 1.0

def compute(Ra):
        consts = [Ra, Pr, S, Pm]
        print(consts)
        #solution = io.fetch_solutions(consts, [0])[0]
        solution = Function(Z)
        d = problem.compute_stability(consts, 0, solution)
        eigs = d.get("eigenvalues", [])
        eigenfunctions = d.get("eigenfunctions", [])
        norm_u = np.array([norm(eigfun.split()[0]) for eigfun in eigenfunctions])
        evals = list(map(complex, eigs))
        evals = np.array([l for l in evals])
        tol = 1.0e-5
        eigsPos = evals[(evals.real>=tol) & (norm_u > tol) & (abs(evals.imag) < 1.0e-12)]
        eigsPos = np.array([l.real for l in eigsPos])
        print(eigsPos)
        print(len(eigsPos))
        outfile = f"output/Ra={Ra}"
        np.save(outfile, eigsPos)
#        plot(Ras[:i+1])
def fix_plot(Ras):
    sol_vec = []
    for i, Ra in enumerate(Ras):
        if i > 0 and i < len(Ras)-2:
            Rasm1 = np.load(f"output/Ra={Ras[i-1]}.npy")
            Ras0 = np.load(f"output/Ra={Ras[i]}.npy")
            Rasp1 = np.load(f"output/Ra={Ras[i+1]}.npy")
#        import ipdb; ipdb.set_trace()
        if (i > 0 and i < len(Ras)-2) and (len(Rasm1) == len(Rasp1)):
            vec = (Rasp1+Rasm1)/2
        else:
            outfile = f"output/Ra={Ra}.npy"
            vec = np.load(outfile)
            
        sol_vec.append(vec[vec<400])
    for xe, ye in zip(Ras, sol_vec):
        plt.scatter([xe] * len(ye), ye, c='b', s=3)
    plot_name = f"eigsPlot_S{int(S)}_Pm{int(Pm)}_B01.png"
    plt.xlabel(r"$\mathrm{Ra}$")
    plt.ylabel(r"$\mathcal{R}(\lambda)$")
    plt.xlim([10**3, 10**5])
    plt.ylim([0, 220])
    plt.savefig(plot_name, dpi=500)

def plot(Ras):
    sol_vec = []
    for Ra in Ras:
        outfile = f"output/Ra={Ra}.npy"
        vec = np.load(outfile)
        sol_vec.append(vec[vec<400])
    for xe, ye in zip(Ras, sol_vec):
        plt.scatter([xe] * len(ye), ye, c='b', s=3)
    plot_name = f"eigsPlot_S{int(S)}_Pm{int(Pm)}_B01.png"
    plt.xlabel(r"$\mathrm{Ra}$")
    plt.ylabel(r"$\mathcal{R}(\lambda)$")
    plt.xlim([10**3, 10**5])
    plt.ylim([0, 220])
    plt.savefig(plot_name, dpi=300)

#compute(100000)
#plot(Ras)

if __name__ == '__main__':
#    pool = Pool(40)
#    pool.map(compute, Ras)
    fix_plot(Ras)
