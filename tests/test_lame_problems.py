import numpy as np
import os
import subprocess
from elastic.meshing import circle, sphere
from elastic import execute, displacement, traction
from errors import check_error, check_interior_error

import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

def circ_from_cart(x, y):
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return r, theta

def sph_from_cart(x, y, z):
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi

def cart_from_circ(r, theta):
    return r * np.cos(theta), r * np.sin(theta)

def cart_from_sph(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


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

def build_trac_bc2d(a, b, p_a, p_b, E, mu):
    def trac_bc(pt):
        r, theta = circ_from_cart(*pt)
        pressure = np.where(r < ((a + b) / 2), -p_a, p_b)
        return cart_from_circ(pressure, theta)
    return trac_bc

def build_trac_bc3d(a, b, p_a, p_b, E, mu):
    def trac_bc(pt):
        r, theta, phi = sph_from_cart(*pt)
        term1 = (p_a * a ** 3 - p_b * b ** 3) / (b ** 3 - a ** 3)
        term2 = ((p_a - p_b) * b ** 3 * a ** 3) / ((b ** 3 - a ** 3) * r ** 3)
        sigmarr = term1 - term2
        sigmarr = np.where(r < ((a + b) / 2.0), sigmarr, -sigmarr)
        return cart_from_sph(sigmarr, theta, phi)
    return trac_bc

def points2d(a, b, n):
    # As a result of the discretization, points on the boundary of the
    # circle that are not vertices in the mesh will be outside the
    # cylinder. Use slightly shifted circle sizes to shift the points
    # inside.
    inside_a = a + 1e-3
    inside_b = b - 1e-3
    nt = int(np.sqrt(n))
    nr = int(np.sqrt(n))
    t_vals = np.linspace(0.0, 2 * np.pi, nt)
    r_vals = np.linspace(inside_a, inside_b, nr)
    r, t = np.meshgrid(r_vals, t_vals)
    x, y = cart_from_circ(r, t)
    return np.array([x.reshape(x.size),y.reshape(y.size)]).T

def points3d(a, b, n):
    # As a result of the discretization, points on the boundary of the
    # circle that are not vertices in the mesh will be outside the
    # cylinder. Use slightly shifted circle sizes to shift the points
    # inside.
    inside_a = a + 1e-1
    inside_b = b - 1e-1
    nt = int(n ** (1.0 / 3.0))
    nphi = int(n ** (1.0 / 3.0))
    nr = int(n ** (1.0 / 3.0))
    r_vals = np.linspace(inside_a, inside_b, nr)
    t_vals = np.linspace(0.0, 2 * np.pi, nt)
    p_vals = np.linspace(0.0, np.pi, nphi)
    r, t, p = np.meshgrid(r_vals, t_vals, p_vals)
    x, y, z = cart_from_sph(r, t, p)
    return np.array([x.reshape(x.size), y.reshape(y.size), z.reshape(z.size)]).T

build_disp_bc = dict()
build_disp_bc[2] = build_disp_bc2d
build_disp_bc[3] = build_disp_bc3d


build_trac_bc = dict()
build_trac_bc[2] = build_trac_bc2d
build_trac_bc[3] = build_trac_bc3d

ball_mesh = dict()
ball_mesh[2] = circle
ball_mesh[3] = sphere

points = dict()
points[2] = points2d
points[3] = points3d

def bc_builder(type_fnc, bc_fnc):
    def f(pts):
        return type_fnc(pts, [bc_fnc(pts[0, :]), bc_fnc(pts[1, :])])
    return f

def lame(dim, bc_types):
    a = 0.8
    b = 1.9
    p_a = 15e6
    p_b = -10e6
    E = 10000.0
    mu = 0.25
    G = E / (2 * (1 + mu))
    refine = dict()
    refine[2] = 6
    refine[3] = 3
    solver_tol = 1e-5

    disp_bc = build_disp_bc[dim](a, b, p_a, p_b, E, mu)
    trac_bc = build_trac_bc[dim](a, b, p_a, p_b, E, mu)
    bc_funcs = dict()
    bc_funcs[traction] = trac_bc
    bc_funcs[displacement] = disp_bc

    es = []
    es.extend(ball_mesh[dim](
        [0] * dim, a, refine[dim],
        bc_builder(bc_types['inner'], bc_funcs[bc_types['inner']]),
        True
    ))
    es.extend(ball_mesh[dim](
        [0] * dim, b, refine[dim],
        bc_builder(bc_types['outer'], bc_funcs[bc_types['outer']]),
        False
    ))
    params = dict(
        shear_modulus = G,
        poisson_ratio = mu,
        solver_tol = solver_tol
    )
    result = execute(dim, es, params)

    import matplotlib.pyplot as plt
    f = result.meshes['continuous'].facets
    xs = f[:, :, 0].reshape(f.shape[0] * f.shape[1])
    ys = f[:, :, 1].reshape(f.shape[0] * f.shape[1])
    fx = result.soln[('continuous', 'traction')][0]
    fy = result.soln[('continuous', 'traction')][1]
    ex, ey = trac_bc([xs, ys])
    # plt.quiver(xs, ys, fx, fy)
    plt.figure()
    plt.plot(xs, fx, 'b')
    plt.plot(xs, fy, 'r')
    plt.figure()
    plt.plot(xs, ex, 'b')
    plt.plot(xs, ey, 'r')
    plt.figure()
    plt.plot(xs, fx - ex, 'b')
    plt.plot(xs, fy - ey, 'r')
    plt.show()

    f = result.meshes['continuous'].facets
    xs = f[:, :, 0].reshape(f.shape[0] * f.shape[1])
    ys = f[:, :, 1].reshape(f.shape[0] * f.shape[1])
    fx = result.soln[('continuous', 'displacement')][0]
    fy = result.soln[('continuous', 'displacement')][1]
    ex, ey = disp_bc([xs, ys])
    # plt.quiver(xs, ys, fx, fy)
    plt.figure()
    plt.plot(xs, fx, 'b')
    plt.plot(xs, fy, 'r')
    plt.figure()
    plt.plot(xs, ex, 'b')
    plt.plot(xs, ey, 'r')
    plt.figure()
    plt.plot(xs, fx - ex, 'b')
    plt.plot(xs, fy - ey, 'r')
    plt.show()

    check_error(result, 'continuous', 'traction', trac_bc, 5e-2)
    check_error(result, 'continuous', 'displacement', disp_bc, 5e-2)

    pts = points[dim](a, b, 1000)
    disp_interior = result.interior_displacement(pts)
    check_interior_error(pts, disp_interior, disp_bc, 3e-2)


def test_disp_disp2d():
    lame(2, dict(inner = displacement, outer = displacement))

def test_trac_trac2d():
    lame(2, dict(inner = traction, outer = traction))

def test_trac_disp2d():
    lame(2, dict(inner = traction, outer = displacement))

def test_disp_trac2d():
    lame(2, dict(inner = displacement, outer = traction))
#
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

if __name__ == '__main__':
    test_disp_disp2d()
