import h5py
import numpy as np
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

def stress(pt):
    b2 = (3.0 / 4.0) * (P / c)
    d4 = -2 * b2 / (c ** 2)
    sxx = -1.5 * (P / (c ** 3)) * pt[0] * pt[1]
    sxy = -0.75 * (P / c) * (1 - (pt[1] ** 2 / c ** 2))
    return [sxx, sxy, sxy, 0]

def trac_bc(pt):
    S = stress(pt)
    if pt[0] < 0.001:
        return S[0], S[2]
    elif pt[0] > 0.99 * L:
        return -S[0], -S[2]
    else:
        return 'bad traction location'


def create_problem():
    refine = 6

    es = []
    es.extend(line([[0, c], [0, -c]], refine, "displacement", disp_bc))
    es.append(Element([[0, -c], [L, -c]], [[0, 0], [0, 0]], "traction", refine))
    es.extend(line([[L, -c], [L, c]], refine, "displacement", disp_bc))
    es.append(Element([[L, c], [0, c]], [[0, 0], [0, 0]], "traction", refine))
    params = dict(
        shear_modulus = G,
        poisson_ratio = nu,
        solver_tol = 5e-5,
        singular_steps = 5,
        obs_order = 2,
        sinh_order = 6,
        dense = False
    )
    return es, params

def test_beam_bend():
    es, params = create_problem()
    result = execute(2, es, params)
    # check_error(result, 'traction', 'displacement', disp_bc, 2e-3)
    # check_error(result, 'displacement', 'traction', trac_bc, 3e-2)

    x, y = np.meshgrid(np.linspace(0, L, 20), np.linspace(-c, c, 20))
    pts = np.array([x.flatten(), y.flatten()]).T
    disp_interior = result.interior_displacement(pts)
    # check_interior_error(pts, disp_interior, disp_bc, 1e-3)

    normals_x = np.array([np.ones(pts.shape[0]), np.zeros(pts.shape[0])]).T
    normals_y = np.array([np.zeros(pts.shape[0]), np.ones(pts.shape[0])]).T
    tx_interior = result.interior_traction(pts, normals_x)
    ty_interior = result.interior_traction(pts, normals_y)
    stress_interior = [tx_interior[0], tx_interior[1], ty_interior[0], ty_interior[1]]
    import matplotlib.pyplot as plt
    plt.quiver(pts[:, 0], pts[:, 1], tx_interior[0], tx_interior[1])
    plt.xlim([-0.1, 1.1])
    plt.ylim([-1.1, 1.1])
    plt.show()
    # check_interior_error(pts, stress_interior, stress, 5e-3)


    f = result.input.meshes['traction'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    s = result.soln[('traction', 'displacement')]
    plt.plot(xs, s[0])
    plt.plot(xs, s[1])
    plt.plot(xs, disp_bc([xs, ys])[0])
    plt.plot(xs, disp_bc([xs, ys])[1])
    plt.show()

    f = result.input.meshes['displacement'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    s = result.soln[('displacement', 'traction')]
    # plt.quiver(xs, ys, s[0], s[1])
    plt.plot(ys, s[0])
    plt.plot(ys, s[1])
    plt.plot(ys, [trac_bc(pt)[0] for pt in zip(xs, ys)])
    plt.plot(ys, [trac_bc(pt)[1] for pt in zip(xs, ys)])
    plt.show()

if __name__ == "__main__":
    test_beam_bend()
