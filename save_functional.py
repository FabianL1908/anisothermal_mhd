from __future__ import absolute_import, print_function

import sys
import os
import gc
import argparse

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
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
rc('text', usetex=True)

from utils import get_branches, get_colors

RB = __import__("rayleigh-benard")

# Branch

#branchids = [44]
#branchids = [34,54,63,64,89]

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--branchids", nargs='+', type=int, default=[-1])
args, _ = parser.parse_known_args()
branchids = args.branchids

#branchids = [0, 38, 2, 42, 166, 161]

comm = COMM_WORLD

# Construct mono-3d problem
problem = RB.RayleighBenardProblem()
mesh = problem.mesh(comm=comm)
Z = problem.function_space(mesh)
functionals = problem.functionals()


if branchids == [-1]:
    branchids = get_branches()
    branchids = [it for el in branchids for item in branchids[el] for it in item]

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
    fixed_params = problem.target_parameter_values()
    knownparams = io.known_parameters(fixed={"Pr": fixed_params["Pr"][0],
                                             "S": fixed_params["S"][0],
                                             "Pm": fixed_params["Pm"][0]}, branchid=branchid)
    knownparams_Ra = np.array([l[0] for l in knownparams])
    Nu = np.array([])
    NT = np.array([])
    NB = np.array([])
    for param in knownparams:
        print(param)
        print("Computing functional for parameters %s, branchid = %d" %
              (str(param[0]), branchid), flush=True)
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

def plot_diagram():
    colors = get_colors()
#    fig = plt.figure()
#    grid = plt.GridSpec(8, 8, hspace=2, wspace=2)
#    fig_u = fig.add_subplot(grid[1:4, 1:4])
#    fig_T = fig.add_subplot(grid[1:4, 4:])
#    fig_B = fig.add_subplot(grid[4:, 4:])
#    figures = [fig_u, fig_T, fig_B]
    for idx, dgrm_type in enumerate(["u", "T", "B"]):
        branchid_dict = get_branches()
        for b_key in branchid_dict:
            color = next(colors)
            for outer_list in branchid_dict[b_key]:
                for branchid in outer_list:
                    with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r') as f:
                        data = list(csv.reader(f, delimiter=","))
                    data = np.array(data)
                    data = data.astype(np.float32)
                    data = data.T
                    plt.plot(data[0], data[1], color=color)
            plt.xlabel(r"$\mathrm{Ra}$")
            plt.ylabel(functionals[idx][2], rotation=0, labelpad=15)
            plt.ylim(bottom=0)
            plt.tight_layout()
            plt.savefig(f'diagram_{dgrm_type}.png', dpi=400)

    colors = get_colors()
    fig = plt.figure()
    grid = plt.GridSpec(4, 4, hspace=2, wspace=2)
    fig_u = fig.add_subplot(grid[:2, :2])
    fig_T = fig.add_subplot(grid[:2, 2:])
    fig_B = fig.add_subplot(grid[2:4, 1:3])
    figures = [fig_u, fig_T, fig_B]
    for idx, dgrm_type in enumerate(["u", "T", "B"]):
        branchid_dict = get_branches()
        for b_key in branchid_dict:
            color = next(colors)
            for outer_list in branchid_dict[b_key]:
                for branchid in outer_list:
                    with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r') as f:
                        data = list(csv.reader(f, delimiter=","))
                    data = np.array(data)
                    data = data.astype(np.float32)
                    data = data.T
                    figures[idx].plot(data[0], data[1], color=color)
            figures[idx].set_xlabel(r"$\mathrm{Ra}$")
            figures[idx].set_ylabel(functionals[idx][2], rotation=0, labelpad=15)
            figures[idx].set_ylim(bottom=0)
#    plt.tight_layout()
        figures[idx].set_xticks(np.linspace(), np.max(data[0]), 5))
    plt.savefig(f'diagram_uTB.png', dpi=400, bbox_inches='tight')
            

plot_diagram()
#shutil.make_archive("/home/boulle/Documents/diagram_data", 'zip', "diagram_data")
