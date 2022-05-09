from __future__ import absolute_import, print_function

import sys
import os
import gc
import csv
import argparse
from functools import partial
from multiprocessing import Pool

from firedrake import *
from firedrake.mg.utils import get_level
import numpy as np
from defcon import *
from defcon import backend, StabilityTask
from defcon.cli.common import fetch_bifurcation_problem
from petsc4py import PETSc

#RB = __import__("rayleigh-benard")
RB = __import__("linear_eigenvalue")

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--branchids", nargs='+', type=int, default=[1])
args, _ = parser.parse_known_args()
branchids = args.branchids

targetparams = np.linspace(0,10**5,401)
minp = 48500
maxp = 55000
targetparams = np.array([t for t in targetparams if (t>=minp and t<=maxp)])

path = "CSV/%d"

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

def get_known_params(branchid):
    # Parameters
    fixed_params = problem.target_parameter_values()
    knownparams = io.known_parameters(fixed={"Pr": fixed_params["Pr"],
                                             "S": fixed_params["S"],
                                             "Pm": fixed_params["Pm"]}, branchid=branchid)
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

    for i in range(0, 10):
        reals = np.array([my_data[param][i][0] for param in params])
        imags = np.array([my_data[param][i][1] for param in params])
        params = np.array(params)
         
        np.savetxt(f"{path_stab}/{branchid}_real_{i}.csv", np.vstack((params, reals)).T, delimiter=",")            
        np.savetxt(f"{path_stab}/{branchid}_imag_{i}.csv", np.vstack((params, imags)).T, delimiter=",")            
    
# Branch
#branchids = [44]
# branchids = [12]
# branchids = [34]
# branchids = [54]
# branchids = [63]
# branchids = [64]
#stab_computation(branchids)
if __name__ == "__main__":
    pool = Pool(4)
    for branchid in branchids:
        knownparams = get_known_params(branchid)
        pool.map(partial(stab_computation, branchid), knownparams)
        create_stability_figures(branchid)
