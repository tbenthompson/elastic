import numpy as np
from elastic import execute, line, displacement, traction
from errors import check_error, check_interior_error
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

# The analytic solution is from the Timoshenko and Goodier book pg. 35
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

def disp_bc(pt, normal):
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

def stress(pt, empty):
    b2 = (3.0 / 4.0) * (P / c)
    d4 = -2 * b2 / (c ** 2)
    sxx = -(P * pt[0] * pt[1]) / I
    sxy = -(P / (2 * I)) * (c ** 2 - pt[1] ** 2)
    return [sxx, sxy, sxy, 0]

def trac_bc(pt, normal):
    S = stress(pt, [])
    tx = S[0] * normal[0] + S[1] * normal[1]
    ty = S[2] * normal[0] + S[3] * normal[1]
    return tx, ty

def create_problem(refine):
    es = []
    es.extend(line([[0, c], [0, -c]], refine,
        lambda pts: displacement(pts, [disp_bc(pts[0, :], []), disp_bc(pts[1, :], [])])
    ))
    es.extend(line([[0, -c], [L, -c]], refine,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    ))
    es.extend(line([[L, -c], [L, c]], refine,
        lambda pts: displacement(pts, [disp_bc(pts[0, :], []), disp_bc(pts[1, :], [])])
    ))
    es.extend(line([[L, c], [0, c]], refine,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    ))
    params = dict(
        shear_modulus = G,
        poisson_ratio = nu,
        length_scale = 0.05,
        solver_tol = 1e-4,
        singular_steps = 8,
        obs_far_order = 5,
        obs_near_order = 5,
        src_far_order = 4,
        sinh_order = 7,
        dense = True,
        check_condition_number = False
    )
    return es, params

def check_soln(result):
    check_error(result, 'continuous', 'displacement', disp_bc, 2e-3)
    check_error(result, 'continuous', 'traction', trac_bc, 1e-2)

    x, y = np.meshgrid(np.linspace(0, L, 20), np.linspace(-c, c, 20))
    pts = np.array([x.flatten(), y.flatten()]).T
    normals_x = np.array([np.ones(pts.shape[0]), np.zeros(pts.shape[0])]).T
    normals_y = np.array([np.zeros(pts.shape[0]), np.ones(pts.shape[0])]).T

    disp_interior = result.interior_displacement(pts)
    check_interior_error(pts, normals_x, disp_interior, disp_bc, 1e-3)

    tx_interior = result.interior_traction(pts, normals_x)
    ty_interior = result.interior_traction(pts, normals_y)
    stress_interior = [tx_interior[0], tx_interior[1], ty_interior[0], ty_interior[1]]
    check_interior_error(pts, normals_x, stress_interior, stress, 4e-3)

def test_beam_bend():
    es, params = create_problem(6)
    result = execute(2, es, params)
    check_soln(result)

if __name__ == "__main__":
    test_beam_bend()
