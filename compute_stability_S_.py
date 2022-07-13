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
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
rc('text', usetex=True)
rc('lines', linewidth=0.7)
rc('lines', markersize=3)
from matplotlib.offsetbox import *
import matplotlib.image as mpimg
from PIL import Image
import warnings
from scipy import ndimage

from firedrake import *
from firedrake.mg.utils import get_level
import numpy as np
from defcon import *
from defcon import backend, StabilityTask
from defcon.cli.common import fetch_bifurcation_problem
from petsc4py import PETSc

from utils import get_branches, get_colors, get_rot_degree_dict, get_image_dict, get_xybox

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

num_eigs = 10

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
                                             "Ra": fixed_params["Ra"][0],
                                             "Pm": fixed_params["Pm"][0]}, branchid=branchid)
    knownparams_init = np.array([l[2] for l in knownparams])
    # if branchid == 284:
    #     knownparams = np.array([t for t in knownparams if (t>=36100)])
#        params = knownparams_init[np.isclose(targetparams[:,None],knownparams_init).any(0)]
    params = knownparams_init
    knownparams = [p for p in knownparams if p[2] in params]
    #print(knownparams)
    return knownparams[::3]

def stab_computation(branchid, param):
#    for param in [knownparams[0]]:
    if os.path.exists(path%(branchid)+f"/{int(param[2])}.csv"):
        print(f"Stability for branchid {branchid} and param {param} was already computed")
    else:
        print("Computing stability for parameters %s, branchid = %d" % (str(param[2]), branchid),flush=True)
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
        np.savetxt(path%(branchid)+"/%.f.csv"%int(param[2]), x, delimiter=",")

def create_pictures():
    for branchid in branchids:
        params = get_known_params(branchid)
        # only save num_plot plots
        num_plots = 20
        len_params = len(params)
        params_idxs = np.linspace(0, len(params)-1, num_plots)
        params_idxs = [int(np.floor(x)) for x in params_idxs]
        if not os.path.isdir("paraview"):
            os.makedirs("paraview")
        for idx in params_idxs:
            solution = io.fetch_solutions(params[idx], [branchid])[0]
            if not os.path.isdir(f"paraview/{branchid}"):
                os.makedirs(f"paraview/{branchid}")
            pvd = File(f"paraview/{branchid}/{int(params[idx][2])}.pvd")
            problem.save_pvd(solution, pvd, params)

    print("Run the following command on your local machine:")
    download_path = "laakmann@vanisher:" + os.getcwd()
    pic_branches = []
    for branchid in branchids:
        mypath = os.path.join("paraview", str(branchid))
        png_files = [f for f in os.listdir(mypath) if os.path.join(mypath, f).endswith(".png")]
        if not len(png_files) > 0:
            pic_branches.append(branchid)
    print(f"mkdir paraview; for branchid in {' '.join([str(num) for num in pic_branches])}; do scp -r {download_path}/paraview/$branchid paraview; done; scp {download_path}/create_png.py .; /Applications/ParaView-5.10.1.app/Contents/bin/pvpython create_png.py; rm -rf paraview/*/*.pvd; rm -rf paraview/*/*.vtu; scp -r paraview/* {download_path}/paraview; rm -rf paraview/*", flush=True)
    input("Hit Enter when you are done:")
        
def create_stability_figures(branchid):
    params = get_known_params(branchid)
    params = [param[2] for param in params]
    path_stab = "StabilityFigures"
    if not os.path.isdir(path_stab):
        os.makedirs(path_stab)

    my_data = {}

    for param in params:
            with open(f'CSV/{branchid}/{int(param)}.csv', 'r') as f:
                    data = list(csv.reader(f, delimiter=","))
            my_data[param] = np.array(data).astype(np.float32)

    try:
        for i in range(0, num_eigs):
            reals = np.array([my_data[param][i][0] for param in params])
            imags = np.array([my_data[param][i][1] for param in params])
            params = np.array(params)

            np.savetxt(f"{path_stab}/{branchid}_real_{i}.csv", np.vstack((params, reals)).T, delimiter=",")            
            np.savetxt(f"{path_stab}/{branchid}_imag_{i}.csv", np.vstack((params, imags)).T, delimiter=",")
    except IndexError:
            print(f"Less than {num_eigs} eigenvalues found")

def get_data(path):
    with open(path, 'r') as f:
        data = list(csv.reader(f, delimiter=","))
    data = np.array(data)
    data = data.astype(np.float32)
    data = data.T
    return data

def add_annotationbox(im_path, x, y, rot_degree):
    l = len(x)
    zipped_list = list(zip(x,y))
    sort_key = lambda x: x[0]
    zipped_list = sorted(zipped_list, key=sort_key)
    xxx = np.array([int(im.split("/")[-1].split("_")[0]) for im in im_path])
    xmid = (x[-1] - x[0]) / 2 + x[0]
    midind = np.abs((xxx-xmid)).argmin()
#    indices = (0, midind-1, len(im_path)-1)
#    import ipdb; ipdb.set_trace()
    image_dict = get_image_dict()
    im_branches = [b.split('/')[1]+b.split('_')[1][0] for b in im_path]
    im_branches = list(dict.fromkeys(im_branches))
#    import ipdb; ipdb.set_trace()
    if im_branches[0] in image_dict:
        indices = []
        pos = []
        for im_b in im_branches:
            indices += [int(f) for f in image_dict[im_b][::2]]
            pos += [f for f in image_dict[im_b][1::2]]
    else:
        indices = (0, min(2,int(midind/2)), midind-1, int((len(im_path)-midind)/2+midind-1), len(im_path)-1)
    if len(im_path) <= indices[-1]:
        im_list = [im_path[i] for i in indices if i < len(im_path)-1]
    else:
        im_list = [im_path[i] for i in indices]
    xx = [int(im.split("/")[-1].split("_")[0]) for im in im_list]
#    xy_list = [zipped_list[ind] for i, ind in enumerate(indices) if int(zipped_list[ind][0]) == xx[i]]
    xy_list = [f for f in zipped_list if int(f[0]) in xx]
    xy_list = list(set(xy_list))
    xy_list = sorted(xy_list, key=sort_key)
    ymin = np.min(y); ymax = np.max(y)
    ymid = (ymax+ymin)/2
    xmin = np.min(x); xmax = np.max(x)
    xmid = (xmax+xmin)/2
    ab_list = []
    for i, xy in enumerate(xy_list):
        if im_branches[0] in image_dict:
            xybox = get_xybox(xy, 0.25*(xmid-xmin), 0.5*(ymid-ymin), pos[i])
        else:
            if xy[1] < ymid:
                xybox = (xy[0], xy[1] + 0.5*(ymid-ymin))
            else:
                xybox = (xy[0], xy[1] - 0.5*(ymid-ymin))
        im = mpimg.imread(im_list[i])
        if rot_degree > 0:
            im = ndimage.rotate(im, rot_degree)
        elif rot_degree == -1:
            im = np.fliplr(im)
        elif rot_degree == -2:
            im = np.flipud(im)
        imagebox = OffsetImage(im, zoom=0.025)
        ab = AnnotationBbox(imagebox, xy,
                            frameon=False,
                            pad=0,
                            xybox=xybox,
                            arrowprops=dict(arrowstyle="->, head_width=0.02, head_length=0.02", linewidth=0.5\
      ))
        ab_list.append(ab)
    return ab_list

def smooth_data(xdata, arr):
#    for i in range(1, len(arr)-2):
#        mid = (arr[i+1]-arr[i-1])*0.5+arr[i-1]
#        if arr[i-1] > 1.0e-8 and abs((mid - arr[i])/arr[i-1]) > 0.2:
#            arr[i] = mid
#    for i in range(2, len(arr)-3):
#        mid = (arr[i+1]-arr[i-1])*0.5+arr[i-1]
#        if (abs((arr[i] - arr[i-1]))/(abs(arr[i-1])+0.001) > 0.5) :
#            arr[i] = mid
#    ind = []
#    for i in range(1, len(arr)-2):
#        if (abs((arr[i] - arr[i-1]))/(abs(arr[i-1])+0.001) > 0.1) :
#           ind.append(i)
#    xdata = np.delete(xdata, ind)
#    arr = np.delete(arr, ind)
    return (xdata, arr)

def plot_stability_figures():
    branchids_dict = get_branches()
    rot_degree_dict = get_rot_degree_dict()
    for b_key in branchids_dict:
#        fig = plt.figure(figsize=(9,6))
#        grid = plt.GridSpec(5, 4, hspace=2, wspace=2)
#        fig_u = fig.add_subplot(grid[:2, :2])
#        fig_T = fig.add_subplot(grid[:2, 2:])
#        fig_B = fig.add_subplot(grid[2:4, 1:3])
#        fig_stab_real = fig.add_subplot(grid[4:, :2])
#        fig_stab_imag = fig.add_subplot(grid[4:, 2:])
#        colors = get_colors()
        len_branch = len(branchids_dict[b_key])
        if len_branch == 1:
            fig = plt.figure(figsize=(10,8))
            grid = plt.GridSpec(10, 8, hspace=14, wspace=14)
            fig_u = fig.add_subplot(grid[:4, :4])
            fig_T = fig.add_subplot(grid[:4, 4:])
            fig_B = fig.add_subplot(grid[4:8, 2:6])
            fig_stab_real = fig.add_subplot(grid[8:, :4])
            fig_stab_imag = fig.add_subplot(grid[8:, 4:])
        elif len_branch == 2 or len_branch==3:
            fig = plt.figure(figsize=(12,12))
            grid = plt.GridSpec(12, 8, hspace=14, wspace=14)
            fig_u = fig.add_subplot(grid[:4, :4])
            fig_T = fig.add_subplot(grid[:4, 4:])
            fig_B = fig.add_subplot(grid[4:8, 2:6])
            fig_stab_real = fig.add_subplot(grid[8:10, :4])
            fig_stab_imag = fig.add_subplot(grid[8:10, 4:])            
            fig_stab_real2 = fig.add_subplot(grid[10:12, :4])
            fig_stab_imag2 = fig.add_subplot(grid[10:12, 4:])            
        else:
            continue
                #raise ValueError("More than two plots per Graph are not possible")
        colors = get_colors()
        color = colors[int(b_key)-1]
        for plot_idx, outer_list in enumerate(branchids_dict[b_key]):
            xdata = np.array([])
            yudata = np.array([])
            yTdata = np.array([])
            yBdata = np.array([])
            from collections import defaultdict
            def def_value():
                return np.array([])
            yrealdata = defaultdict(def_value)
            yimagdata = defaultdict(def_value)
            max_num_branch = 100
            for branchid in outer_list:
                print(f"{b_key}: Plotting branch {branchid}")
                data = get_data(f'diagram_u/{branchid}.csv')
                xdata = np.append(xdata, data[0])
                yudata = np.append(yudata, data[1])
                data = get_data(f'diagram_T/{branchid}.csv')
                yTdata = np.append(yTdata, data[1])
                data = get_data(f'diagram_B/{branchid}.csv')
                yBdata = np.append(yBdata, data[1])
                for i in range(0, num_eigs):
                    if i >= max_num_branch:
                        continue
                    try:
                        data = get_data(f'StabilityFigures/{branchid}_real_{i}.csv')
                        yrealdata[i] = np.append(yrealdata[i], data[1])
                        data = get_data(f'StabilityFigures/{branchid}_imag_{i}.csv')
                        yimagdata[i] = np.append(yimagdata[i], data[1])
                    except FileNotFoundError:
                        max_num_branch = i
                        print(f"Less than {num_eigs} eigenvalues found")
            argsort = np.argsort(xdata)
            xdata = xdata[argsort]
            yudata = yudata[argsort]
            yTdata = yTdata[argsort]
            yBdata = yBdata[argsort]
            for key in yrealdata:
                yrealdata[key] = yrealdata[key][argsort]
            for key in yimagdata:
                yimagdata[key] = yimagdata[key][argsort]
            fig_u.plot(xdata, yudata, color=color)
            fig_T.plot(xdata, yTdata, color=color)
            fig_B.plot(xdata, yBdata, color=color)
            color3 = "b"
            for i in range(0, num_eigs):
                try:
#                    with warnings.catch_warnings():
#                        warnings.simplefilter("ignore", category=RuntimeWarning)
                    if len(yrealdata[i]) == 0 or np.mean(yrealdata[i]) > 400:
                        continue
                    xdata_real, smooth_yrealdata = smooth_data(xdata, yrealdata[i]) 
                    xdata_imag, smooth_yimagdata = smooth_data(xdata, yimagdata[i]) 
                    if plot_idx == 0:
#                        if b_key == '1':
#                            import ipdb; ipdb.set_trace()
                        fig_stab_real.scatter(xdata_real, smooth_yrealdata, color=color3, marker='.')
                        fig_stab_imag.scatter(xdata_imag, smooth_yimagdata, color=color3, marker='.')
                    elif plot_idx == 1:
                        fig_stab_real2.scatter(xdata_real, smooth_yrealdata, color=color3, marker='.')
                        fig_stab_imag2.scatter(xdata_imag, smooth_yimagdata, color=color3, marker='.')
                except FileNotFoundError:
                    print(f"Less than {num_eigs} eigenvalues found")
                except ValueError:
                    print(f"Less than {num_eigs} eigenvalues found")
            image_files = [os.listdir(f"paraview/{branch}") for branch in outer_list]
            image_files = [os.path.join(f"paraview/{branch}", f) for branch in outer_list for f in os.listdir(f"paraview/{branch}")]
            image_files = [f for f in image_files if f.endswith(".png")]
            sort_key = lambda x: int(x.split("/")[-1].split("_")[0])
            image_files = sorted(image_files, key=sort_key)
            u_image_files = [f for f in image_files if f.endswith("u.png")]           
            ab_list = add_annotationbox(u_image_files, xdata, yudata, rot_degree_dict[int(b_key)])
            for ab in ab_list:
                fig_u.add_artist(ab)
            T_image_files = [f for f in image_files if f.endswith("T.png")]
            ab_list = add_annotationbox(T_image_files, xdata, yTdata, 0)
            for ab in ab_list:
                fig_T.add_artist(ab)
            B_image_files = [f for f in image_files if f.endswith("B.png")]
            ab_list = add_annotationbox(B_image_files, xdata, yBdata, rot_degree_dict[int(b_key)])
            for ab in ab_list:
                fig_B.add_artist(ab)
        fig_u.set_xlabel(r"$\mathrm{S}$")
#        import matplotlib.ticker as ticker
#        fig_u.set_major_locator(ticker.MultipleLocator(1))
        fig_T.set_xlabel(r"$\mathrm{Ra}$")
        fig_B.set_xlabel(r"$\mathrm{Ra}$")
        fig_u.set_ylabel(problem.functionals()[0][2], rotation=0, labelpad=15)
        fig_T.set_ylabel(problem.functionals()[1][2], rotation=0, labelpad=15)
        fig_B.set_ylabel(problem.functionals()[2][2], rotation=0, labelpad=15)
        xlims = fig_u.get_xlim()
         
        fig_stab_real.set_xlabel(r"$\mathrm{S}$")
        fig_stab_imag.set_xlabel(r"$\mathrm{S}$")
        fig_stab_real.set_ylabel(r"$\mathcal{R}(\lambda)$", rotation=0, labelpad=15)
        fig_stab_imag.set_ylabel(r"$\mathcal{I}(\lambda)$", rotation=0, labelpad=15)
        y0 = fig_stab_real.get_ylim()[0]
        fig_stab_real.set_ylim(bottom=y0-2)
        y0 = fig_stab_imag.get_ylim()[0]
        fig_stab_imag.set_ylim(bottom=y0-2)
        if fig_stab_real.get_ylim()[1] < 10.0:
#            y0 = fig_stab_real.get_ylim()[0]
            fig_stab_real.set_ylim(top=10)
#            fig_stab_real.set_ylim(bottom=y0-2)
        fig_stab_imag.set_ylim(bottom=0)
        fig_stab_imag.set_ylim(top=max(fig_stab_imag.get_ylim()[1], 0.01))
        fig_stab_real.axhline(0, color='black')
        fig_stab_real.set_xlim(xlims)
        fig_stab_imag.set_xlim(xlims)
        if len_branch == 2:
            fig_stab_real2.set_xlabel(r"$\mathrm{S}$")
            fig_stab_imag2.set_xlabel(r"$\mathrm{S}$")
            fig_stab_real2.set_ylabel(r"$\mathcal{R}(\lambda)$", rotation=0, labelpad=15)
            fig_stab_imag2.set_ylabel(r"$\mathcal{I}(\lambda)$", rotation=0, labelpad=15)
            fig_stab_imag2.set_ylim(bottom=0)
            y0 = fig_stab_real2.get_ylim()[0]
            fig_stab_real2.set_ylim(bottom=y0-2)
            y0 = fig_stab_imag2.get_ylim()[0]
            fig_stab_imag2.set_ylim(bottom=y0-2)
            if fig_stab_real2.get_ylim()[1] < 10.0:
#                y0 = fig_stab_real2.get_ylim()[1]
                fig_stab_real2.set_ylim(top=10)
#                fig_stab_real2.set_ylim(bottom=y0-2)
            fig_stab_real2.axhline(0, color='black')
            fig_stab_real2.set_xlim(xlims)
            fig_stab_imag2.set_xlim(xlims)
#        fig_stab_imag.set_ylim(bottom=-0.001)
            real_ylim = list(fig_stab_real.get_ylim())
            imag_ylim = list(fig_stab_imag.get_ylim())
            real_ylim[0] = min(fig_stab_real.get_ylim()[0], fig_stab_real2.get_ylim()[0])
            real_ylim[1] = max(fig_stab_real.get_ylim()[1], fig_stab_real2.get_ylim()[1])
            imag_ylim[0] = min(fig_stab_imag.get_ylim()[0], fig_stab_imag.get_ylim()[0])
            imag_ylim[1] = max(fig_stab_imag.get_ylim()[1], fig_stab_imag.get_ylim()[1])
            fig_stab_real.set_ylim(real_ylim)
            fig_stab_imag.set_ylim(imag_ylim)
            fig_stab_real2.set_ylim(real_ylim)
            fig_stab_imag2.set_ylim(imag_ylim)
            
        plt.savefig(f'StabilityFigures/diagram_branch_{b_key}.png', dpi=800)
    
    
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
    create_pictures()
    plot_stability_figures()
#    for branchid in [188]:
#        knownparams = get_known_params(branchid)
#        pool.map(partial(stab_computation, branchid), knownparams)
#        create_stability_figures(branchid)
#        stab_computation(branchid, knownparams[-1]) 
