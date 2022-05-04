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

RB = __import__("rayleigh-benard")

def stab_computation(branchids, comm=COMM_WORLD):
    #targetparams = np.linspace(0,10**5,201)

    targetparams = np.linspace(0,10**5,401)
    minp = 48500
    maxp = 55000
    targetparams = np.array([t for t in targetparams if (t>=minp and t<=maxp)])
    
    path = "CSV/%d"
    
    #comm = COMM_WORLD
    
    # Construct mono-3d problem
    problem = RB.RayleighBenardProblem()
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
    
    for branchid in branchids:   
        # Parameters
        knownparams = io.known_parameters(fixed={"Pr":1}, branchid=branchid)
        knownparams_init = np.array([l[0] for l in knownparams])
        # if branchid == 284:
        #     knownparams = np.array([t for t in knownparams if (t>=36100)])
        params = knownparams_init[np.isclose(targetparams[:,None],knownparams_init).any(0)]
        knownparams = [p for p in knownparams if p[0] in params]
        print(knownparams)
        for param in knownparams:
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
            
# Branch
#branchids = [44]
# branchids = [12]
# branchids = [34]
# branchids = [54]
# branchids = [63]
# branchids = [64]
branchids = [215]
stab_computation(branchids)


Hi Lisa,

I'm a DPhil student in my final year in Oxford and I'm looking for an accommodation for my last two month in Oxford for July and August. The pictures look really nice and I would be very interested to take either of the two rooms you offer. Can you give me a few more details about the flat and whom I would be living with? And how much you usually pay per month for the additional bills?

I'm looking forward to hearing from you. Best,
Fabian

Thanks for replying so fast. That sounds really great. I'm doing my DPhil in maths, too :-). Would you mind sharing the name of your housemate, maybe I even know him. And I would love to come by to have a look at the room. Just let me know when it's the most convenient for you, I'm quite flexible.

Hi Zhoulin, I was wondering if you are still looking for a housemate. I was just told that I cannot extend my contract with Keble beyond the 30th June and now I'm looking for an accommodation for at least July and August. I've seen some rooms on Facebook but I thought I ask you first in case you are interested. If you don't want to share your flat or feel uncomfortable with it please just say so, that's no problem at all.
 
