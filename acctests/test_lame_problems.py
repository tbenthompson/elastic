import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from tools.input_builder import *
from tools.coordinate_transforms import *

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
        sigmarr = np.where(r > ((a + b) / 2.0), sigmarr, -sigmarr)
        return cart_from_sph(sigmarr, theta, phi)
    return trac_bc

build_trac_bc = dict()
build_trac_bc[2] = build_trac_bc2d
build_trac_bc[3] = build_trac_bc3d

ball_mesh = dict()
ball_mesh[2] = circle
ball_mesh[3] = sphere

def concentric_circle_pts(a, b, nt, nr):
    t_vals = np.linspace(0.0, 2 * np.pi, nt)
    r_vals = np.linspace(a, b, nr)
    r, t = np.meshgrid(r_vals, t_vals)
    return cart_from_circ(r, t)

def plotter(a, b, dim, disp_bc):
    nt = 50
    nr = 20
    pts = list(concentric_circle_pts(a, b, nt, nr))
    if dim == 3:
        pts.append(np.zeros_like(pts[0]))

    soln = disp_bc(*pts)

    x = pts[0]
    y = pts[1]
    ux = soln[0]
    uy = soln[1]
    plt.figure()
    plt.contourf(x, y, ux)
    plt.contour(x, y, ux, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.contourf(x, y, uy)
    plt.contour(x, y, uy, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.quiver(x, y, ux, uy)
    plt.show()

def delete_files(dir, filename):
    root = os.path.splitext(filename)[0]
    for f in os.listdir(dir):
        if root in f:
            os.remove(os.path.join(dir, f))

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

def lame(dim, bc_types):
    a = 0.8
    b = 1.9
    p_a = 10e6
    p_b = -15e6
    E = 80e9
    mu = 0.25
    G = E / (2 * (1 + mu))
    refine = dict()
    refine[2] = 8
    refine[3] = 3
    solver_tol = 1e-10
    dir = 'test_data/auto_gen/'
    filename = 'lame.in'
    input_filepath = dir + filename
    pts_filepath = dir + filename + '_pts'

    disp_bc = build_disp_bc[dim](a, b, p_a, p_b, E, mu)
    # plotter(a, b, dim, disp_bc)
    trac_bc = build_trac_bc[dim](a, b, p_a, p_b, E, mu)
    bc_funcs = dict()
    bc_funcs['traction'] = trac_bc
    bc_funcs['displacement'] = disp_bc

    in_root, file_ext = os.path.splitext(input_filepath)
    displacement_filepath = in_root + '.disp_out'
    traction_filepath = in_root + '.trac_out'
    interior_disp_filepath = in_root + '.disp_out_interior'

    delete_files(dir, filename)

    es = []
    es.extend(ball_mesh[dim]([0] * dim, a, refine[dim], bc_types['inner'],
                     bc_funcs[bc_types['inner']], True))
    es.extend(ball_mesh[dim]([0] * dim, b, refine[dim], bc_types['outer'],
                     bc_funcs[bc_types['outer']], False))

    bem_template(input_filepath, es = es,
                 shear_modulus = G, mu = mu, solver_tol = solver_tol)
    run(input_filepath, dim = dim, stdout_dest = subprocess.PIPE)
    if 'displacement' in bc_types.values():
        check_field(traction_filepath, trac_bc, False, -6)
    if 'traction' in bc_types.values():
        check_field(displacement_filepath, disp_bc, False, 6)

    nt = 20
    nr = 20
    points(a, b, nt, nr, pts_filepath)
    interior_run(input_filepath, pts_filepath)
    check_field(interior_disp_filepath, disp_bc, False, 6)


def test_disp_disp2d():
    lame(2, dict(inner = "displacement", outer = "displacement"))

def test_trac_disp2d():
    lame(2, dict(inner = "traction", outer = "displacement"))

def test_disp_trac2d():
    lame(2, dict(inner = "displacement", outer = "traction"))

def test_trac_trac2d():
    lame(2, dict(inner = "traction", outer = "traction"))

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
