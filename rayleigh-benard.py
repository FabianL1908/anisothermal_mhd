# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 13:25:15 2019

@author: user
"""

# -*- coding: utf-8 -*-
import sys
import os
import re
from math import floor

from petsc4py import PETSc
from firedrake import *
from defcon import *
from matplotlib import pyplot as plt
import numpy
import numpy as np

from defcon.parametertools import parameters_to_string

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

# This hack enforces the boundary condition at (0, 0)


class PointwiseBC(DirichletBC):
    @utils.cached_property
    def nodes(self):
        x = self.function_space().mesh().coordinates.dat.data_ro
        zero = numpy.array([0, 0])
        dists = [numpy.linalg.norm(pt - zero) for pt in x]
        minpt = numpy.argmin(dists)
        if dists[minpt] < 1.0e-10:
            out = numpy.array([minpt], dtype=numpy.int32)
            print("done")
        else:
            out = numpy.array([], dtype=numpy.int32)
        return out


class RayleighBenardProblem(BifurcationProblem):
    def mesh(self, comm):
        #        self.mesh = RectangleMesh(50, 50, 1, 1, quadrilateral=True, comm=comm)
        self.mesh = RectangleMesh(50, 50, 1, 1, quadrilateral=False, diagonal="crossed", comm=comm)
        self._solver = None
        return self.mesh

    def function_space(self, mesh):
        k = 2
        V = VectorFunctionSpace(mesh, "CG", k)
        Q = FunctionSpace(mesh, "CG", k-1)
        #Vel = FiniteElement("N2div", mesh.ufl_cell(), k, variant="integral")
        #V = FunctionSpace(mesh, Vel)
        # Q = FunctionSpace(mesh, "DG", k-1)  # p
        R = FunctionSpace(mesh, "CG", k)  # E
        Wel = FiniteElement("N1div", mesh.ufl_cell(), k, variant="integral")
        W = FunctionSpace(mesh, Wel)
        TT = FunctionSpace(mesh, "CG", k)
        self.Z = MixedFunctionSpace([V, Q, TT, W, R])
        return self.Z

    def parameters(self):
        Ra = Constant(0.0)
        Pr = Constant(0.0)
        S = Constant(0.0)
        Pm = Constant(0.0)
        return [
            (Ra, "Ra", r"$\mathrm{Ra}$"),
            (Pr, "Pr", r"$\mathrm{Pr}$"),
            (S, "S", r"$\mathrm{S}$"),
            (Pm, "Pm", r"$\mathrm{Pm}$"),
        ]

    def residual(self, z, params, w):
        (Ra, Pr, S, Pm) = params
        (u, p, T, B, E) = split(z)
        (v, q, s, C, Ff) = split(w)
        g = Constant((0, 1))
        nn = FacetNormal(self.mesh)
        gamma = Constant(1)
        eps = lambda x: sym(grad(x))

        F = (
              2 * inner(eps(u), eps(v))*dx
            + inner(dot(grad(u), u), v) * dx
            + gamma * inner(div(u), div(v)) * dx
            + S * inner(vcross(B, E), v) * dx
            + S * inner(vcross(B, scross(u, B)), v) * dx
            - inner(p, div(v)) * dx
            - inner(div(u), q) * dx
            + inner(E, Ff) * dx
            + inner(scross(u, B), Ff) * dx
            - 1/Pm * inner(B, vcurl(Ff)) * dx
            + inner(vcurl(E), C) * dx
            + 1/Pm * inner(div(B), div(C)) * dx
            - Ra/Pr * inner(g*T, v) * dx
            + 1/Pr * inner(grad(T), grad(s)) * dx
            + inner(dot(u, grad(T)), s) * dx
#            - inner(dot(grad(T), nn), s)*ds(3)
#            - inner(dot(grad(T), nn), s)*ds(4)
            )

#        F = (
#              inner(grad(u), grad(v))*dx
#            + inner(dot(u,grad(u)),v)*dx
#            - inner(p, div(v))*dx
#            - Ra*Pr*inner(T*g, v)*dx
#            + inner(div(u), q)*dx
#            + inner(grad(T), grad(S))*dx
#            + inner(dot(grad(T), u), S)*dx
#            - inner(dot(grad(T),nn),S)*ds(3)
#            - inner(dot(grad(T),nn),S)*ds(4)
#          )
        return F

    
    def boundary_conditions(self, Z, params):
        u0, p0, T0, B0, E0 = self.equilibrium_solution()        
        bcs = [
                DirichletBC(Z.sub(0), u0, (1, 2, 3, 4)), # u = 0 at the Boundary
                DirichletBC(Z.sub(2), Constant(1.0), (3)), # T = 1 at y = 0
                DirichletBC(Z.sub(2), Constant(0.0), (4)), # T = 0 at y = pi
                DirichletBC(Z.sub(3), B0, "on_boundary"),
                DirichletBC(Z.sub(4), E0, "on_boundary"),
                ]
        return bcs

    def functionals(self):
        def u_L2(z, params):
            (u, p, T, B, E) = split(z)
            j = assemble(inner(u, u)*dx)
            return j
        
        def T_L2(z, params):
            (u, p, T, B, E) = split(z)
            j = assemble(inner(T, T)*dx)
            return j

        def B_L2(z, params):
            (u, p, T, B, E) = split(z)
            j = assemble(inner(B, B)*dx)
            return j

        return [(u_L2, "u_sqL2", r"$\|u\|^2$"),
                (T_L2, "T_sqL2", r"$\|T\|^2$"),
                (B_L2, "B_sqL2", r"$\|B\|^2$")]

    def get_initial_indices(self, params):
#        params = [1.0e5, 1.0, 100.0, 1.0]
        mypath = "./output/" + parameters_to_string(self.parameters(), params)
#        print(mypath)
        myfiles = os.listdir(mypath)
        indices = [re.findall('solution-\d*', f) for f in myfiles]
        indices = [int(ind[0].split("-")[1]) for ind in indices if len(ind)>0]
        indices.sort()
        return indices
    
    def number_initial_guesses(self, params):
        init_guesses = self.get_initial_indices(params)
        print(f"Number of initial guesses: {len(init_guesses)}")
        print(f"Initial indices: {init_guesses}")
        return len(init_guesses)

    def initial_guess(self, Z, params, n):
        consts = [params[0], params[1], params[2], params[3]]
        print(f"consts = {consts}")
        print(type(consts[0]))
        params = self.parameters()
#        params[3] = (Constant(1.0), params[3][1], params[3][2])
        a = self.io("output")
        a.setup(params, self.functionals(), Z)
#        Listi = range(self.number_initial_guesses(params))
        indices = self.get_initial_indices(consts)
        print(indices)
        sols = a.fetch_solutions(consts, indices)

        return sols[n]

    def equilibrium_solution(self):
        x = SpatialCoordinate(self.mesh)
        u0 = Constant((0.0, 0.0), domain=self.mesh)
        p0 = Constant(0.0, domain=self.mesh)
        T0 = 1.0 - x[1]
        B0 = Constant((0.0, 1.0), domain=self.mesh)
        E0 = Constant(0.0, domain=self.mesh)

        return (u0, p0, T0, B0, E0)

    def trivial_solutions(self, Z, params, freeindex):
        z = Function(self.Z)
        x = SpatialCoordinate(self.mesh)

        u0, p0, T0, B0, E0 = self.equilibrium_solution()
        z.sub(0).interpolate(u0)
        z.sub(1).interpolate(p0)
        z.sub(2).interpolate(T0)
        z.sub(3).interpolate(B0)
        z.sub(4).interpolate(E0)
        return [z]

    def number_solutions(self, params):
        (Ra, Pr, S, Pm) = params
        #if Ra < 1700:
        #    return 1
        #if Ra < 1720:
        #    return 3
        return float("inf")

    def squared_norm(self, z, w, params):
        (zu, zp, zT, zB, zE) = split(z)
        (wu, wp, wT, wB, wE) = split(w)
        diffu = zu - wu
        diffp = zp - wp
        diffT = zT - wT
        diffB = zB - wB
        diffE = zE - wE
        
        diff = (
              inner(diffu, diffu)*dx
            + inner(grad(diffu), grad(diffu))*dx
            + inner(diffT, diffT)*dx
            + inner(diffB, diffB)*dx
            )
        return diff

    def save_pvd(self, z, pvd, params):
        (u, p, T, B, E) = z.split()
        u.rename("Velocity", "Velocity")
        p.rename("Pressure", "Pressure")
        T.rename("Temperature", "Temperature")
        B.rename("MagneticField", "MagneticField")
        E.rename("ElectricFieldf", "ElectricFieldf")
        pvd.write(u, p, T, B, E)
    
    def solver_parameters(self, params, task, **kwargs):
        linesearch = "basic"
        damping = 1.0

        lu = {
             "mat_type": "aij",
             "snes_max_it": 30,
             "snes_type": "newtonls",
             "snes_linesearch_type": linesearch,
             "snes_stol": 0.0,
             "snes_atol": 1.0e-7,
             "snes_rtol": 0.0,
#             "snes_divergence_tolerance": 1.0e4,
             "snes_monitor": None,
             "snes_converged_reason": None,
             "snes_linesearch_monitor": None,
             "ksp_type": "preonly",
             "ksp_monitor_true_residual": None,
             "ksp_max_it": 10,
             "pc_type": "lu",
             "pc_factor_mat_solver_type": "mumps",
             }
        return lu

    def solver(self, problem, params, solver_params, prefix="", **kwargs):
        gc.collect()
        # Check for failed PCs before
        if self._solver is not None:
            if self._solver.snes.ksp.getConvergedReason() == PETSc.KSP.ConvergedReason.DIVERGED_PCSETUP_FAILED:
                print("Resetting solver")
                self._solver = None

        if self._solver is None:
            nullspace = MixedVectorSpaceBasis(self.Z, [self.Z.sub(0), VectorSpaceBasis(constant=True), self.Z.sub(2), self.Z.sub(3), self.Z.sub(4)])
            self._solver = NonlinearVariationalSolver(problem, nullspace = nullspace, options_prefix=prefix,
                                 solver_parameters=solver_params)
        return self._solver

    def parameter_values(self):
        values = {#"Ra": linspace(10**5, 10**3, 40),
                  "Ra": [10**5],
                  "Pr": [1.0],
                  "S": [1.0],
                  "Pm": [1.0],
        }
        return values

if __name__ == "__main__":
    if False:
        dc = DeflatedContinuation(problem=RayleighBenardProblem(), teamsize=1, verbose=True, disable_deflation=True)
        values = {#"Ra": linspace(10**5, 10**3, 40),
                         "Ra": [10**5],
                         "Pr": [1.0],
                         "S": [1.0],
                         "Pm": linspace(1.0, 10.0, 20),
            }
        dc.run(values=values, freeparam="Pm")
    if False:
        dc = DeflatedContinuation(problem=RayleighBenardProblem(), teamsize=1, verbose=True, disable_deflation=True)
        values = {#"Ra": linspace(10**5, 10**3, 40),
                         "Ra": [10**5],
                         "Pr": [1.0],
                         "S": linspace(1.0, 1000.0, 20),
                         "Pm": [10.0],
            }
        dc.run(values=values, freeparam="S")
    if True:
        dc = DeflatedContinuation(problem=RayleighBenardProblem(), teamsize=1, verbose=True)
        values = {"Ra": linspace(10**5, 10**3, 50),
                          "Pr": [1.0],
                          "S": [1000.0],
                          "Pm": [10.0],
            }
        dc.run(values=values, freeparam="Ra")
