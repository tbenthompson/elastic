import subprocess
import os
import numpy as np
from elastic import line, crack, free_slip, execute
from errors import check_error

import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.INFO)

def build_slip_bc(a, stress_drop, G, nu):
    def slip_bc(pt, normal):
        # Classic Griffith crack solution. Normally, the solution is given
        # as a displacement for each crack face. I multiply by 2 to get the
        # slip across the crack.
        slip = -2 * (1 - nu) * (stress_drop / G) * np.sqrt(a ** 2 - pt[0] ** 2)
        return slip, np.zeros_like(pt[0])
    return slip_bc

def create_problem(a, stress_drop, G, nu, bc_type):
    refine = 6
    es = []
    es.extend(line([[-a, 0], [a, 0]], refine, lambda pts: bc_type(pts, stress_drop)))
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
    problem = execute(2, es, params)
    check_error(problem, 'discontinuous', 'slip', slip_bc, 4e-2)

def test_crack():
    griffith_soln(lambda pts, bc: crack(pts, [[bc, 0], [bc, 0]]))

def test_free_slip():
    griffith_soln(lambda pts, bc: free_slip(pts, [0, 0], [[bc, bc]]))

if __name__ == '__main__':
    test_crack()
