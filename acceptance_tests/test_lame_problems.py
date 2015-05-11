import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from elastic.mesh_gen import *
from elastic.solver import Controller
from coordinate_transforms import *

def build_disp_bc2d(a, b, p_a, p_b, E, mu):
    def disp_bc(pt):
        r, theta = circ_from_cart(*pt)
        ur = ((1 + mu) * a ** 2 * b ** 2) / (E * (b ** 2 - a ** 2)) *\
            (((p_a - p_b) / r) +
             ((1 - 2 * mu) * ((p_a * a ** 2 - p_b * b ** 2) / (a ** 2 * b ** 2)) * r))
        return cart_from_circ(ur, theta)
    return disp_bc

def build_disp_bc3d(a, b, p_a, p_b, E, mu):
    def disp_bc(pt):
        r, theta, phi = sph_from_cart(*pt)
        ur = (1.0 / (2 * E * (b ** 3 - a ** 3) * r ** 2)) *\
            (2 * (p_a * a ** 3 - p_b * b ** 3) * (1 - 2 * mu) * r ** 3 +
             (p_a - p_b) * (1 + mu) * b ** 3 * a ** 3)
        return cart_from_sph(ur, theta, phi)
    return disp_bc

build_disp_bc = dict()
build_disp_bc[2] = build_disp_bc2d
build_disp_bc[3] = build_disp_bc3d

def build_trac_bc2d(a, b, p_a, p_b, E, mu):
    def trac_bc(pt):
        r, theta = circ_from_cart(*pt)
        pressure = np.where(r < ((a + b) / 2), p_a, -p_b)
        return cart_from_circ(pressure, theta)
    return trac_bc

def build_trac_bc3d(a, b, p_a, p_b, E, mu):
    def trac_bc(pt):
        r, theta, phi = sph_from_cart(*pt)
        term1 = (p_a * a ** 3 - p_b * b ** 3) / (b ** 3 - a ** 3)
        term2 = ((p_a - p_b) * b ** 3 * a ** 3) / ((b ** 3 - a ** 3) * r ** 3)
        sigmarr = term1 - term2
        sigmarr = np.where(r > ((a + b) / 2.0), -sigmarr, sigmarr)
        return cart_from_sph(sigmarr, theta, phi)
    return trac_bc

build_trac_bc = dict()
build_trac_bc[2] = build_trac_bc2d
build_trac_bc[3] = build_trac_bc3d

ball_mesh = dict()
ball_mesh[2] = circle
ball_mesh[3] = sphere

def points(a, b, nt, nr, out_filepath):
    # As a result of the discretization, points on the boundary of the
    # circle that are not vertices in the mesh will be outside the
    # cylinder. Use slightly shifted circle sizes to shift the points
    # inside.
    inside_a = a + 1e-3
    inside_b = b - 1e-3
    x, y = concentric_circle_pts(inside_a, inside_b, nt, nr)
    pts = zip(x.flatten(), y.flatten())
    points_template(out_filepath, pts)

def l2_error(dim, mesh, exact_fnc, est_vec):
    l2_diff = np.zeros(dim)
    l2_exact = np.zeros(dim)
    for f_idx in range(mesh.facets.shape[0]):
        for v_idx in range(dim):
            global_v_idx = f_idx * dim + v_idx
            v = mesh.facets[f_idx, v_idx, :]
            exact = exact_fnc(v)
            for d_est in range(dim):
                est = est_vec[d_est][global_v_idx]
                l2_diff[d_est] += (exact[d_est] - est) ** 2
                l2_exact += exact[d_est] ** 2
    error = np.sqrt(l2_diff / l2_exact)
    return error

def check_error(problem, mesh_name, field_name, exact_fnc, limit):
    if problem.input.meshes[mesh_name].facets.shape[0] > 0:
        error = l2_error(problem.dim,
            problem.input.meshes[mesh_name],
            exact_fnc,
            problem.soln[(mesh_name, field_name)]
        )
        print('(mesh = ' + mesh_name + ', field = ' + field_name + ') error: '\
             + str(error))
        assert(np.all(error < limit))

def lame(dim, bc_types):
    a = 0.8
    b = 1.9
    p_a = 15e6
    p_b = -10e6
    E = 80e9
    mu = 0.25
    G = E / (2 * (1 + mu))
    refine = dict()
    refine[2] = 9
    refine[3] = 3
    solver_tol = 1e-5

    disp_bc = build_disp_bc[dim](a, b, p_a, p_b, E, mu)
    trac_bc = build_trac_bc[dim](a, b, p_a, p_b, E, mu)
    bc_funcs = dict()
    bc_funcs['traction'] = trac_bc
    bc_funcs['displacement'] = disp_bc

    es = []
    es.extend(ball_mesh[dim]([0] * dim, a, refine[dim], bc_types['inner'],
                     bc_funcs[bc_types['inner']], True))
    es.extend(ball_mesh[dim]([0] * dim, b, refine[dim], bc_types['outer'],
                     bc_funcs[bc_types['outer']], False))
    params = dict(
        shear_modulus = G,
        poisson_ratio = mu,
        solver_tol = solver_tol
    )
    problem = Controller(dim, es, params)

    # m = problem.input.meshes['traction']
    # s = problem.soln[('traction', 'displacement')]
    # verts = m.facets[:,:,:].reshape((m.facets.shape[0] * dim, dim))
    # plt.figure()
    # plt.plot(verts[:,0], s[1], 'k.')
    # plt.figure()
    # plt.plot(verts[:,0], [disp_bc(v)[1] for v in verts], 'b.')
    # plt.show()

    check_error(problem, 'displacement', 'traction', trac_bc, 4e-2)
    check_error(problem, 'traction', 'displacement', disp_bc, 4e-2)

    # nt = 20
    # nr = 20
    # points(a, b, nt, nr, pts_filepath)
    # interior_run(input_filepath, pts_filepath)
    # check_field(interior_disp_filepath, disp_bc, False, 6)


def test_disp_disp2d():
    lame(2, dict(inner = "displacement", outer = "displacement"))

def test_trac_trac2d():
    lame(2, dict(inner = "traction", outer = "traction"))

def test_trac_disp2d():
    lame(2, dict(inner = "traction", outer = "displacement"))

def test_disp_trac2d():
    lame(2, dict(inner = "displacement", outer = "traction"))

# def test_disp_disp3d():
#     lame(3, dict(inner = "displacement", outer = "displacement"))
#
# def test_trac_disp3d():
#     lame(3, dict(inner = "traction", outer = "displacement"))
#
# def test_disp_trac3d():
#     lame(3, dict(inner = "displacement", outer = "traction"))
#
# def test_trac_trac3d():
#     lame(3, dict(inner = "traction", outer = "traction"))
