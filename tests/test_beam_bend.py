import h5py
import numpy as np
import subprocess
from elastic import execute, line, displacement, traction
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

def stress(pt):
    b2 = (3.0 / 4.0) * (P / c)
    d4 = -2 * b2 / (c ** 2)
    sxx = -1.5 * (P / (c ** 3)) * pt[0] * pt[1]
    sxy = -0.75 * (P / c) * (1 - (pt[1] ** 2 / c ** 2))
    return [sxx, sxy, sxy, 0]

def trac_bc(pt):
    S = stress(pt)
    if pt[0] == 0:
        return S[0], S[2]
    elif pt[0] == L:
        return -S[0], -S[2]


def create_problem():
    refine = 6

    es = []
    es.extend(line([[0, c], [0, -c]], refine,
        lambda pts: displacement(pts, [disp_bc(pts[0, :]), disp_bc(pts[1, :])])
    ))
    es.extend(line([[0, -c], [L, -c]], refine,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    ))
    es.extend(line([[L, -c], [L, c]], refine,
        lambda pts: displacement(pts, [disp_bc(pts[0, :]), disp_bc(pts[1, :])])
    ))
    es.extend(line([[L, c], [0, c]], refine,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    ))
    params = dict(
        shear_modulus = G,
        poisson_ratio = nu,
        solver_tol = 1e-6,
        singular_steps = 8,
        obs_order = 4,
        sinh_order = 7,
        dense = True
    )
    return es, params

def test_beam_bend():
    es, params = create_problem()
    result = execute(2, es, params)
    check_error(result, 'continuous', 'displacement', disp_bc, 2e-3)

    x, y = np.meshgrid(np.linspace(0, L, 20), np.linspace(-c, c, 20))
    pts = np.array([x.flatten(), y.flatten()]).T
    disp_interior = result.interior_displacement(pts)
    check_interior_error(pts, disp_interior, disp_bc, 1e-3)

    normals_x = np.array([np.ones(pts.shape[0]), np.zeros(pts.shape[0])]).T
    normals_y = np.array([np.zeros(pts.shape[0]), np.ones(pts.shape[0])]).T
    tx_interior = result.interior_traction(pts, normals_x)
    ty_interior = result.interior_traction(pts, normals_y)
    stress_interior = [tx_interior[0], tx_interior[1], ty_interior[0], ty_interior[1]]
    check_interior_error(pts, stress_interior, stress, 4e-3)

if __name__ == "__main__":
    test_beam_bend()
