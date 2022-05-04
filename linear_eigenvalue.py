# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 13:25:15 2019

@author: user
"""

# -*- coding: utf-8 -*-
import sys
from   math import floor

from petsc4py import PETSc
from slepc4py import SLEPc

from firedrake import *
from defcon import *

import numpy
import numpy as np

# Define two-dimensional versions of cross and curl operators
def scross(x, y):
    return x[0]*y[1] - x[1]*y[0]


def vcross(x, y):
    return as_vector([x[1]*y, -x[0]*y])


def scurl(x):
    return x[1].dx(0) - x[0].dx(1)


def vcurl(x):
    return as_vector([x.dx(1), -x.dx(0)])


def acurl(x):
    return as_vector([
                     x[2].dx(1),
                     -x[2].dx(0),
                     x[1].dx(0) - x[0].dx(1)
                     ])

import importlib
RB = importlib.import_module("rayleigh-benard")

class EVRayleighBenardProblem(RB.RayleighBenardProblem):
    def boundary_conditions(self, Z, params):
        bcs = [
                DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (1, 2, 3, 4)), # u = 0 at the boundary
                DirichletBC(Z.sub(2), Constant(0.0), (3, 4)), # T = 0 at y = 0 and y = pi 3/4
                DirichletBC(Z.sub(3), Constant((0.0, 0.0)), (1, 2, 3, 4)),
                DirichletBC(Z.sub(4), Constant(0.0), (1, 2, 3, 4)),
                ]
        return bcs

    def number_initial_guesses(self, params):
        return 1

    def initial_guess(self, Z, params, n):
        x = SpatialCoordinate(Z.mesh())
        z = Function(Z)
        return z
    
    def compute_stability(self, params, branchid, z, hint=None):                
        (Ra, Pr, S, Pm) = params
        
        trial = TrialFunction(self.Z)
        (u, p, T, B, E) = split(trial)
        test = TestFunction(self.Z)
        (v, q, s, C, Ff) = split(test)
        
        g = Constant((0,1))
        nn  = FacetNormal(self.mesh)
        gamma = Constant(0.0) #self.constants()["gamma"]
        u0, p0, T0, B0, E0 = self.equilibrium_solution()
        x = SpatialCoordinate(self.mesh)
        eps = lambda x: sym(grad(x))

        stabform = (
            - 2 * inner(eps(u), eps(v))*dx
            + Ra*Pr*inner(T*g, v)*dx
            - 1/Pr * inner(grad(T), grad(s))*dx
            + inner(-dot(grad(T0), u), s)*dx
            - gamma * inner(div(u), div(v)) * dx
            - S * inner(vcross(B0, E), v) * dx
            - S * inner(vcross(B0, scross(u, B0)), v) * dx
            + inner(p, div(v)) * dx
            + inner(div(u), q) * dx
            - inner(E, Ff) * dx
            - inner(scross(u, B0), Ff) * dx
            + 1/Pm * inner(B, vcurl(Ff)) * dx
            - inner(vcurl(E), C) * dx
            - 1/Pm * inner(div(B), div(C)) * dx
        )

        my_z = Function(self.Z)
        my_z.split()[0].interpolate(u0)
        my_z.split()[1].interpolate(p0)
        my_z.split()[2].interpolate(T0)
        my_z.split()[3].interpolate(B0)
        my_z.split()[4].interpolate(E0)
        Fsol = -self.residual(my_z, params, test)
        stabform = derivative(Fsol, my_z, trial)
        massform = inner(u, v)*dx + inner(T, s)*dx + inner(B, C)*dx #+ inner(E, Ff)*dx
#        stabform -= Ra*Pr*inner(split(trial)[2]*g, v)*dx 
#        massform = -Pr*inner(T*g, v)*dx + inner(B,C)*dx

        stabbcs = self.boundary_conditions(self.Z, params)
        M = assemble(massform, bcs=stabbcs, mat_type="aij")

        # There must be a better way of doing this
        from firedrake.preconditioners.patch import bcdofs
        for bc in stabbcs:
            M.M.handle.zeroRowsColumns(bcdofs(bc), diag=0.0)
        stabmass = M

        comm = self.Z.mesh().comm
        
        A = assemble(stabform, bcs=stabbcs, mat_type="aij")

        # Create the SLEPc eigensolver
        eps = SLEPc.EPS().create(comm=comm)
        eps.setOperators(A.M.handle, stabmass.M.handle)
        eps.setProblemType(eps.ProblemType.GNHEP)
        eps.setFromOptions()
        
        eps.solve()
        eigenvalues = []
        eigenfunctions = []
        eigenfunction = Function(self.Z, name="Eigenfunction")

        for i in range(eps.getConverged()):
            lmbda = eps.getEigenvalue(i)
            eigenvalues.append(lmbda)
            with eigenfunction.dat.vec_wo as x:
                eps.getEigenvector(i, x)
            eigenfunctions.append(eigenfunction.copy(deepcopy=True))

        tol = 1.0e-10
        if any(abs(omega.imag) > tol for omega in eigenvalues):
            is_stable = False
        else:
            is_stable = True

        # Shift and normalize the eigenfunctions
        #x = SpatialCoordinate(self.mesh)
        #shift = Function(self.Z)
        #u0, p0, T0, B0, E0 = self.equilibrium_solution()
        #shift.sub(0).interpolate(u0)
        #shift.sub(1).interpolate(p0)
        #shift.sub(2).interpolate(T0)
        #shift.sub(3).interpolate(B0)
        #shift.sub(4).interpolate(E0)
        #for eigenfunction in eigenfunctions:
        #    eigenfunction.assign(eigenfunction+shift)

        d = {"stable": is_stable,
             "eigenvalues": eigenvalues,
             "eigenfunctions": eigenfunctions,
             "hint": eigenfunctions}

        return d

    def solver_parameters(self, params, task, **kwargs):
        linesearch = "basic"
        damping = 1.0

        lu = {
             "mat_type": "aij",
             "snes_max_it": 100,
             "snes_type": "newtonls",
             "snes_linesearch_type": linesearch,
             "snes_stol": 0.0,
             "snes_atol": 1.0e-13,
             "snes_rtol": 0.0,
             "snes_divergence_tolerance": -1,
             "snes_monitor": None,
             "snes_converged_reason": None,
             "snes_linesearch_monitor": None,
             "ksp_type": "preonly",
             "ksp_monitor_true_residual": None,
             "ksp_max_it": 10,
             "pc_type": "lu",
             "pc_factor_mat_solver_type": "mumps",
#             "eps_monitor_all": None,
             "eps_converged_reason": None,
             "eps_type": "krylovschur",
             "eps_nev": 15,
             "eps_max_it": 30,
             "eps_tol": 1e-8,
             "eps_target": 500,
             "epswhich": "smallest_magnitude",
             "st_type": "sinvert",
             "st_ksp_type": "preonly",
             "st_pc_type": "lu",
             "st_pc_factor_mat_solver_type": "mumps",
             "st_mat_mumps_icntl_14": 5000,
             "st_ksp_max_it": 10,
             "ds_parallel": "synchronized"
             }
        return lu

    def solver(self, problem, params, solver_params, prefix="", **kwargs):
        self._solver = NonlinearVariationalSolver(problem, options_prefix=prefix, solver_parameters=solver_params)
        return self._solver
