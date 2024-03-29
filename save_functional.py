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
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
plt.rcParams['ytick.right'] = True
labelpad = 30

from utils import get_colors, get_linestyles

RB = __import__("rayleigh-benard")

# Branch

#branchids = [44]
#branchids = [34,54,63,64,89]

def get_branches():
    branch_dict = {}
    with open('branches3.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            branch_dict[row[0]] = []
    with open('branches3.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            int_row = [int(r) for r in row[1:]]
            branch_dict[row[0]].append(int_row)
    return branch_dict

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--branchids", nargs='+', type=int, default=[-1])
parser.add_argument("--mode", choices=["Ra", "S"], type=str, required=True)
args, _ = parser.parse_known_args()
branchids = args.branchids
mode = args.mode

#branchids = [0, 38, 2, 42, 166, 161]

comm = COMM_WORLD

if mode == "Ra":
    IDX = 0
elif mode == "S":
    IDX = 2

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

def save_functional():
    for branchid in branchids:
        # Parameters
        fixed_params = problem.target_parameter_values()
        if mode == "Ra":
            knownparams = io.known_parameters(fixed={"Pr": fixed_params["Pr"][0],
                                                     "S": fixed_params["S"][0],
                                                     "Pm": fixed_params["Pm"][0]}, branchid=branchid)
        elif mode == "S":
            knownparams = io.known_parameters(fixed={"Pr": fixed_params["Pr"][0],
                                                     "Ra": fixed_params["Ra"][0],
                                                     "Pm": fixed_params["Pm"][0]}, branchid=branchid)[::3]
        knownparams_S = np.array([l[IDX] for l in knownparams])
        Nu = np.array([])
        NT = np.array([])
        NB = np.array([])
        for param in knownparams:
            print(param)
            print("Computing functional for parameters %s, branchid = %d" %
                  (str(param[IDX]), branchid), flush=True)
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
        knownparams_S = knownparams_S.reshape((len(knownparams_S), 1))
        Nu = Nu.reshape((len(Nu), 1))
        NT = NT.reshape((len(NT), 1))
        NB = NB.reshape((len(NB), 1))
        np.savetxt("diagram_u/%d.csv"%branchid, np.hstack((knownparams_S, Nu)), delimiter=",")
        np.savetxt("diagram_T/%d.csv"%branchid, np.hstack((knownparams_S, NT)), delimiter=",")
        np.savetxt("diagram_B/%d.csv"%branchid, np.hstack((knownparams_S, NB)), delimiter=",")

def get_data(dgrm_type, branchid):
    with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r') as f:
        data = list(csv.reader(f, delimiter=","))
    data = np.array(data)
    return data

def write_data(dgrm_type, branchid, left, targetp, avg_val):
    line = f"{targetp},{avg_val}"
    if left:
        with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)
    else:
        with open(f'diagram_{dgrm_type}/{branchid}.csv', 'a') as f:
            f.write(line)

def join_plots():
    try:
        with open('join_plots.csv', 'r') as f:
            data = list(csv.reader(f, delimiter=","))
    except:
        return

    for dat in data:
        branchid_1, branchid_2, left, scale = dat
        left = bool(int(left))
        scale = float(scale)
        for dgrm_type in ["u", "T", "B"]:
            data1 = get_data(dgrm_type, branchid_1)
            data2 = get_data(dgrm_type, branchid_2) 
            endp1, val1 = data1[0] if left else data1[-1]
            endp2, val2 = data2[0] if left else data2[-1]
            if  endp1 != endp2:
                endp1 = (float(endp1)+float(endp2))/2
                endp2 = endp1
            targetp = float(endp1) - scale*np.abs(float(data1[1][0]) - float(data1[0][0])) if left else float(endp1) + scale*np.abs(float(data1[-1][0]) - float(data1[-2][0]))
            avg_val = (float(val1) + float(val2))/2
            write_data(dgrm_type, branchid_1, left, targetp, avg_val)
            write_data(dgrm_type, branchid_2, left, targetp, avg_val)
        

def plot_diagram():
#    fig = plt.figure()
#    grid = plt.GridSpec(8, 8, hspace=2, wspace=2)
#    fig_u = fig.add_subplot(grid[1:4, 1:4])
#    fig_T = fig.add_subplot(grid[1:4, 4:])
#    fig_B = fig.add_subplot(grid[4:, 4:])
#    figures = [fig_u, fig_T, fig_B]
    for idx, dgrm_type in enumerate(["u", "T", "B"]):
#        import ipdb; ipdb.set_trace()
        plt.figure()
        colors = get_colors()
        linestyles = get_linestyles()
        branchid_dict = get_branches()
        for b_key in branchid_dict:
            color = colors[int(b_key)-1]
            linestyle = linestyles[int(b_key)-1]
            for outer_list in branchid_dict[b_key]:
                full_data = np.array([]).reshape((0,2))
                for branchid in outer_list:
                    print(branchid)
                    with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r') as f:
                        data = list(csv.reader(f, delimiter=","))
                    data = np.array(data)
                    data = data.astype(np.float32)
                    full_data = np.vstack((full_data, data))
                full_data = full_data[full_data[:, 0].argsort()]
                full_data = full_data.T
                plt.plot(full_data[0], full_data[1], color=color, label=f"{b_key}", linestyle=linestyle)
        xlabel_str = r"$\mathrm{" + mode + "}$"
        plt.xlabel(xlabel_str)
        plt.ylabel(functionals[idx][2], rotation=0, labelpad=labelpad)
        if dgrm_type == "u":
            plt.ylim(bottom=0)
        if dgrm_type == "B":
            plt.ylim(bottom=1)
        right = 10**3 if mode == "S" else 10**5
        plt.xlim(right=right)
        plt.xlim(left=0)
#        plt.tight_layout()
        handles, labels = plt.gca().get_legend_handles_labels()
        if dgrm_type == "B":
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.0, 1.0))
 #            plt.legend(bbox_to_anchor=(1, 1.0))
        if dgrm_type == "T":
            plt.hlines(y=1/3, xmin=0, xmax=10**5, color='black', linewidth=0.6, linestyle="--")
        plt.savefig(f'diagram_{dgrm_type}.png', dpi=400, bbox_inches="tight")

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
            color = colors[int(b_key)-1]
            for outer_list in branchid_dict[b_key]:
                for branchid in outer_list:
                    with open(f'diagram_{dgrm_type}/{branchid}.csv', 'r') as f:
                        data = list(csv.reader(f, delimiter=","))
                    data = np.array(data)
                    data = data.astype(np.float32)
                    data = data.T
                figures[idx].plot(data[0], data[1], color=color)
            figures[idx].set_xlabel(xlabel_str)
            figures[idx].set_ylabel(functionals[idx][2], rotation=0, labelpad=labelpad)
            figures[idx].set_ylim(bottom=0)
            right = 10**3 if mode == "S" else 10**5
            figures[idx].set_xlim(right=right)
            figures[idx].set_xlim(left=0)
#    plt.tight_layout()
#        figures[idx].set_xticks(np.linspace(), np.max(data[0]), 5))
    plt.savefig(f'diagram_uTB.png', dpi=400, bbox_inches='tight')
            
#save_functional()
#join_plots()
plot_diagram()
#shutil.make_archive("/home/boulle/Documents/diagram_data", 'zip', "diagram_data")
