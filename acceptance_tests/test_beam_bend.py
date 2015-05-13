import h5py
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from elastic.mesh_gen import line
from elastic.solver import execute
from elastic.input_builder import Element
from errors import check_error, check_interior_error

def G_from_E_nu(E, nu):
    return E / (2 * (1 + nu))

L = 1.0
c = 1.0
I = (2.0 / 3.0) * c ** 3

P = -40e6
E = 80e9
nu = 0.25
G = G_from_E_nu(E, nu)

# transformation from plane stress to plane strain
nu_fic = nu / (1 - nu)
E_fic = E / (1 - nu ** 2)
G_fic = G_from_E_nu(E_fic, nu_fic)

def disp_bc(pt):
    x = pt[0]
    y = pt[1]
    ux = (-P * x ** 2 * y) / (2 * E_fic * I) \
          - (nu_fic * P * y ** 3) / (6 * E_fic * I) \
          + (P * y ** 3) / (6 * I * G_fic) \
          + (P * L ** 2 * y) / (2 * E_fic * I) \
          - (P * c ** 2 * y) / (2 * I * G_fic)

    uy = (nu_fic * P * x * y ** 2) / (2 * E_fic * I) \
        + (P * x ** 3) / (6 * E_fic * I) \
        - (P * L ** 2 * x) / (2 * E_fic * I) \
        + (P * L ** 3) / (3 * E_fic * I)

    return ux, uy

def create_problem():
    refine = 5

    es = []
    es.extend(line([[0, c], [0, -c]], refine, "displacement", disp_bc))
    es.append(Element([[0, -c], [L, -c]], [[0, 0], [0, 0]], "traction", refine))
    es.extend(line([[L, -c], [L, c]], refine, "displacement", disp_bc))
    es.append(Element([[L, c], [0, c]], [[0, 0], [0, 0]], "traction", refine))
    params = dict(
        shear_modulus = G,
        poisson_ratio = nu,
    )
    return es, params

def test_beam_bend():
    es, params = create_problem()
    result = execute(2, es, params)
    check_error(result, 'traction', 'displacement', disp_bc, 2e-3)

    x, y = np.meshgrid(np.linspace(0, L, 20), np.linspace(-c, c, 20))
    pts = np.array([x.flatten(), y.flatten()]).T
    disp_interior = result.interior_displacement(pts)
    check_interior_error(pts, disp_interior, disp_bc, 1e-3)

if __name__ == "__main__":
    plotter()
