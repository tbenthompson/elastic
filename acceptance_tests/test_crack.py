import subprocess
import os
import numpy as np
from elastic.solver import Controller
from elastic.input_builder import Element
from errors import check_error

def build_slip_bc(a, stress_drop, G, nu):
    def slip_bc(pt):
        # Classic Griffith crack solution. Normally, the solution is given
        # as a displacement for each crack face. I multiply by 2 to get the
        # slip across the crack.
        slip = 2 * (1 - nu) * (stress_drop / G) * np.sqrt(a ** 2 - pt[0] ** 2)
        return slip, np.zeros_like(pt[0])
    return slip_bc

def create_problem(a, stress_drop, G, nu, bc_type):
    refine = 6
    bc = [[stress_drop, 0], [stress_drop, 0]]
    es = []
    es.append(Element([[-a, 0], [-a/1.3, 0]], bc, bc_type, refine))
    es.append(Element([[-a/1.3, 0], [a/1.3, 0]], bc, bc_type, refine))
    es.append(Element([[a/1.3, 0], [a, 0]], bc, bc_type, refine))
    params = dict(
        shear_modulus = G,
        poisson_ratio = nu
    )

    return es, params

def griffith_soln(bc_type):
    a = 0.5
    G = 40e9
    nu = 0.29
    stress_drop = 5.5e6

    slip_bc = build_slip_bc(a, stress_drop, G, nu)
    es, params = create_problem(a, stress_drop, G, nu, bc_type)
    problem = Controller(2, es, params)
    check_error(problem, 'crack_traction', 'slip', slip_bc, 4e-2)

def test_crack():
    griffith_soln("crack_traction")

def test_free_slip():
    griffith_soln("free_slip_traction")
