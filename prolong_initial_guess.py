from firedrake import *
from defcon import *
import numpy as np
from alfi import *

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
branchid = 6
fixed_params = problem.target_parameter_values()
knownparams = io.known_parameters(fixed={"Pr": fixed_params["Pr"][0],
                                             "S": fixed_params["S"][0],
                                             "Pm": fixed_params["Pm"][0]}, branchid=branchid)

sols = io.fetch_solutions(knownparams[0], [6])
sol = sols[0]
mesh1 = problem.mesh
mesh2 = MeshHierarchy(mesh1, 1)[-1]
Z2 = problem.function_space(mesh2)
qtransfer = NullTransfer()
dgtransfer = DGInjection()
Etransfer = DGInjection()

transfers = {
             Z.sub(1).ufl_element(): (prolong, restrict, qtransfer.inject),
             Z.sub(4).ufl_element(): (prolong, restrict, Etransfer.inject),
             VectorElement("DG", mesh.ufl_cell(), 1): (dgtransfer.prolong, 
restrict, dgtransfer.inject),
            }
tm = TransferManager(native_transfers=transfers)

sol_prol = Function(Z2)
tm.prolong(sol, sol_prol)
io2 = problem.io("output2")
io2.setup(problem.parameters(), problem.functionals(), Z2)
import ipdb; ipdb.set_trace()
#functionals = io2.fetch_functionals(knownparams, branchid)
io2.save_solution(sol_prol, [], knownparams[0], branchid)
pvd = File("initial_guess/solution/guess.pvd")
problem.save_pvd(sol, pvd, knownparams[0])
problem.save_pvd(sol_prol, pvd, knownparams[0])
