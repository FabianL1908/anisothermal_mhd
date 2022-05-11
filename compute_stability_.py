from __future__ import absolute_import, print_function

import sys
import os
import gc
import csv
import argparse
from functools import partial
from multiprocessing import Pool
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
rc('text', usetex=True)

from firedrake import *
from firedrake.mg.utils import get_level
import numpy as np
from defcon import *
from defcon import backend, StabilityTask
from defcon.cli.common import fetch_bifurcation_problem
from petsc4py import PETSc

from utils import get_branches, get_colors

#RB = __import__("rayleigh-benard")
RB = __import__("linear_eigenvalue")

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--branchids", nargs='+', type=int, default=[-1])
args, _ = parser.parse_known_args()
branchids = args.branchids

targetparams = np.linspace(0, 10**5, 401)
minp = 48500
maxp = 55000
targetparams = np.array([t for t in targetparams if (t >= minp and t <= maxp)])

path = "CSV/%d"

comm = COMM_WORLD

# Construct mono-3d problem
problem = RB.EVRayleighBenardProblem()
mesh = problem.mesh(comm=comm)
Z = problem.function_space(mesh)

if branchids == [-1]:
    branchids = get_branches()
    branchids = [item for el in branchids for item in el]


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

def get_known_params(branchid):
    # Parameters
    fixed_params = problem.target_parameter_values()
    knownparams = io.known_parameters(fixed={"Pr": fixed_params["Pr"][0],
                                             "S": fixed_params["S"][0],
                                             "Pm": fixed_params["Pm"][0]}, branchid=branchid)
    knownparams_init = np.array([l[0] for l in knownparams])
    # if branchid == 284:
    #     knownparams = np.array([t for t in knownparams if (t>=36100)])
#        params = knownparams_init[np.isclose(targetparams[:,None],knownparams_init).any(0)]
    params = knownparams_init
    knownparams = [p for p in knownparams if p[0] in params]
    print(knownparams)
    return knownparams    

def stab_computation(branchid, param):
#    for param in [knownparams[0]]:
    if os.path.exists(path%(branchid)+f"/{int(param[0])}.csv"):
        print(f"Stability for branchid {branchid} and param {param} was already computed")
    else:
        print("Computing stability for parameters %s, branchid = %d" % (str(param[0]), branchid),flush=True)

        consts = param
        #try:

        solution = io.fetch_solutions(consts, [branchid])[0]
        d = problem.compute_stability(consts, branchid, solution)
        evals = list(map(complex, d["eigenvalues"]))
        RpointsMu = np.array([l.real for l in evals])
        RpointsMu = RpointsMu.reshape(len(RpointsMu),1)
        IpointsMu = np.array([l.imag for l in evals])
        IpointsMu = IpointsMu.reshape(len(IpointsMu),1)
        x = np.hstack((RpointsMu,IpointsMu))

        # Sort x by largest real part
        x = np.flipud(x[np.lexsort(np.fliplr(x).T)])

        # Save the eigenvalues
        if not os.path.isdir(path%(branchid)):
            os.makedirs(path%(branchid))
        np.savetxt(path%(branchid)+"/%.f.csv"%param[0], x, delimiter=",")

def create_stability_figures(branchid):
    params = get_known_params(branchid)
    params = [param[0] for param in params]
    path_stab = "StabilityFigures"
    if not os.path.isdir(path_stab):
        os.makedirs(path_stab)

    my_data = {}

    for param in params:
            with open(f'CSV/{branchid}/{int(param)}.csv', 'r') as f:
                    data = list(csv.reader(f, delimiter=","))
            my_data[param] = np.array(data).astype(np.float32)

    try:
        for i in range(0, 10):
            reals = np.array([my_data[param][i][0] for param in params])
            imags = np.array([my_data[param][i][1] for param in params])
            params = np.array(params)

            np.savetxt(f"{path_stab}/{branchid}_real_{i}.csv", np.vstack((params, reals)).T, delimiter=",")            
            np.savetxt(f"{path_stab}/{branchid}_imag_{i}.csv", np.vstack((params, imags)).T, delimiter=",")
    except IndexError:
            print("Less than 10 eigenvalues found")

def get_data(path):
    with open(path, 'r') as f:
        data = list(csv.reader(f, delimiter=","))
    data = np.array(data)
    data = data.astype(np.float32)
    data = data.T
    return data

def plot_stability_figures():
    branchids_list = get_branches()
    for branchid_l in branchids_list:
        fig = plt.figure()
        grid = plt.GridSpec(5, 4, hspace=2, wspace=2)
        fig_u = fig.add_subplot(grid[:2, :2])
        fig_T = fig.add_subplot(grid[:2, 2:])
        fig_B = fig.add_subplot(grid[2:4, 1:3])
        fig_stab_real = fig.add_subplot(grid[4:, :2])
        fig_stab_imag = fig.add_subplot(grid[4:, 2:])
        colors = get_colors()
        color = next(colors)
        for branchid in branchid_l:
            data = get_data(f'diagram_u/{branchid}.csv')
            fig_u.plot(data[0], data[1], color=color)
            data = get_data(f'diagram_T/{branchid}.csv')
            fig_T.plot(data[0], data[1], color=color)
            data = get_data(f'diagram_B/{branchid}.csv')
            fig_B.plot(data[0], data[1], color=color)
            colors2 = get_colors()
            for i in range(0, 10):
                color2 = next(colors2)
                try:
                    data = get_data(f'StabilityFigures/{branchid}_real_{i}.csv')
                    fig_stab_real.plot(data[0], data[1], color=color2)
                    data = get_data(f'StabilityFigures/{branchid}_imag_{i}.csv')
                    fig_stab_imag.plot(data[0], data[1], color=color2)
                except FileNotFoundError:
                    print("Less than 10 eigenvalues found")
        fig_u.set_xlabel(r"$\mathrm{Ra}$")
        fig_T.set_xlabel(r"$\mathrm{Ra}$")
        fig_B.set_xlabel(r"$\mathrm{Ra}$")
        fig_stab_real.set_xlabel(r"$\mathrm{Ra}$")
        fig_stab_imag.set_xlabel(r"$\mathrm{Ra}$")
        fig_u.set_ylabel(problem.functionals()[0][2])
        fig_T.set_ylabel(problem.functionals()[1][2])
        fig_B.set_ylabel(problem.functionals()[2][2])
        fig_stab_real.set_ylabel(r"$\mathcal{R}(\lambda)$")
        fig_stab_imag.set_ylabel(r"$\mathcal{I}(\lambda)$")
        plt.savefig(f'diagram_branch_{branchid_l[0]}.png', dpi=400)
    
    
# Branch
#branchids = [44]
# branchids = [12]
# branchids = [34]
# branchids = [54]
# branchids = [63]
# branchids = [64]
#stab_computation(branchids)
if __name__ == "__main__":
    pool = Pool(40)
    print(branchids)
    for branchid in branchids:
        knownparams = get_known_params(branchid)
        pool.map(partial(stab_computation, branchid), knownparams)
        create_stability_figures(branchid)
    plot_stability_figures()
