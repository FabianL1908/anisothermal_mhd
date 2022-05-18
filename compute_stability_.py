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
from matplotlib.offsetbox import *
import matplotlib.image as mpimg
from PIL import Image


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
    branchids = [it for el in branchids for item in branchids[el] for it in item]


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
    #print(knownparams)
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

def create_pictures():
    for branchid in branchids:
        params = get_known_params(branchid)
        # only save 10 plots
        len_params = len(params)
        params_idxs = np.linspace(0, len(params)-1, 10)
        params_idxs = [int(np.floor(x)) for x in params_idxs]
        if not os.path.isdir("paraview"):
            os.makedirs("paraview")
        for idx in params_idxs:
            solution = io.fetch_solutions(params[idx], [branchid])[0]
            if not os.path.isdir(f"paraview/{branchid}"):
                os.makedirs(f"paraview/{branchid}")
            pvd = File(f"paraview/{branchid}/{int(params[idx][0])}.pvd")
            problem.save_pvd(solution, pvd, params)

    print("The the following command on your local machine:")
    download_path = "laakmann@wolverine:" + os.getcwd()
    print(f"scp -r {download_path}/paraview .; scp {download_path}/create_png.py .; /Applications/ParaView-5.10.1.app/Contents/bin/pvpython create_png.py; rm -rf paraview/*/*.pvd; rm -rf paraview/*/*.vtu; scp -r paraview/* {download_path}/paraview; rm -rf paraview/*")
    input("Hit Enter when you are done:")
        
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

def add_annotationbox(im_path, x, y):
    l = len(x)

    zipped_list = list(zip(x,y))
    sort_key = lambda x: x[0]
    zipped_list = sorted(zipped_list, key=sort_key)
    indices = (0, 4, 9)
    im_list = [im_path[i] for i in indices]
    xx = [int(im.split("/")[-1].split("_")[0]) for im in im_list]
#    xy_list = [zipped_list[ind] for i, ind in enumerate(indices) if int(zipped_list[ind][0]) == xx[i]]
    xy_list = [f for f in zipped_list if f[0] in xx]

    ymin = np.min(y); ymax = np.max(y)
    ymid = (ymax+ymin)/2
    ab_list = []
    for i, xy in enumerate(xy_list):
        if xy[1] < ymid:
            xybox = (xy[0], xy[1] + 0.5*(ymid-ymin))
        else:
            xybox = (xy[0], xy[1] - 0.5*(ymid-ymin))
        im = mpimg.imread(im_list[i])
        imagebox = OffsetImage(im, zoom=0.02)
        ab = AnnotationBbox(imagebox, xy,
                            frameon=False,
                            pad=0,
                            xybox=xybox,
                            arrowprops=dict(arrowstyle="->, head_width=0.02, head_length=0.02", linewidth=0.5\
      ))
        ab_list.append(ab)
    return ab_list

def plot_stability_figures():
    branchids_dict = get_branches()
    for b_key in branchids_dict:
        fig = plt.figure()
        grid = plt.GridSpec(5, 4, hspace=2, wspace=2)
        fig_u = fig.add_subplot(grid[:2, :2])
        fig_T = fig.add_subplot(grid[:2, 2:])
        fig_B = fig.add_subplot(grid[2:4, 1:3])
        fig_stab_real = fig.add_subplot(grid[4:, :2])
        fig_stab_imag = fig.add_subplot(grid[4:, 2:])
        colors = get_colors()
        for outer_list in branchids_dict[b_key]:
            xdata = np.array([])
            yudata = np.array([])
            yTdata = np.array([])
            yBdata = np.array([])
            color = next(colors)
            for branchid in outer_list:
                data = get_data(f'diagram_u/{branchid}.csv')
                xdata = np.append(xdata, data[0])
                yudata = np.append(yudata, data[1])
                fig_u.plot(data[0], data[1], color=color)
                data = get_data(f'diagram_T/{branchid}.csv')
                yTdata = np.append(yTdata, data[1])
                fig_T.plot(data[0], data[1], color=color)
                data = get_data(f'diagram_B/{branchid}.csv')
                yBdata = np.append(yBdata, data[1])
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
#            image_files = [os.listdir(f"paraview/{branch}") for outer_list in branchids_dict[b_key] for branch in outer_list]
            image_files = [os.listdir(f"paraview/{branch}") for branch in outer_list]
#            image_files = [os.path.join(f"paraview/{branch}", f) for outer_list in branchids_dict[b_key] for branch in outer_list for f in os.listdir(f"paraview/{branch}")]
            image_files = [os.path.join(f"paraview/{branch}", f) for branch in outer_list for f in os.listdir(f"paraview/{branch}")]
            image_files = [f for f in image_files if f.endswith(".png")]
            sort_key = lambda x: int(x.split("/")[-1].split("_")[0])
            image_files = sorted(image_files, key=sort_key)
            u_image_files = [f for f in image_files if f.endswith("u.png")]
            ab_list = add_annotationbox(u_image_files, xdata, yudata)
            for ab in ab_list:
                fig_u.add_artist(ab)
            T_image_files = [f for f in image_files if f.endswith("T.png")]
            ab_list = add_annotationbox(T_image_files, xdata, yTdata)
            for ab in ab_list:
                fig_T.add_artist(ab)
            B_image_files = [f for f in image_files if f.endswith("B.png")]
            ab_list = add_annotationbox(B_image_files, xdata, yBdata)
            for ab in ab_list:
                fig_B.add_artist(ab)
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
        
        plt.savefig(f'diagram_branch_{b_key}.png', dpi=400)
    
    
# Branch
#branchids = [44]
# branchids = [12]
# branchids = [34]
# branchids = [54]
# branchids = [63]
# branchids = [64]
#stab_computation(branchids)
if __name__ == "__main__":
#    pool = Pool(40)
    print(branchids)
#    for branchid in branchids:
#        knownparams = get_known_params(branchid)
#        pool.map(partial(stab_computation, branchid), knownparams)
#        create_stability_figures(branchid)
#    create_pictures()
    plot_stability_figures()
