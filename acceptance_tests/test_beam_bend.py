import h5py
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from elastic.mesh_gen import line
from elastic.solver import Controller
from elastic.input_builder import Element
from errors import check_error

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

def upper_lower_trac_bc(x, y):
    return np.zeros_like(x), np.zeros_like(x)

def plotter():
    nx = 20
    ny = 20
    x_vals = np.linspace(0.0, L, nx)
    y_vals = np.linspace(-c, c, ny)
    y_vals = [-c, c]
    x, y = np.meshgrid(x_vals, y_vals)

    ux, uy = disp_bc(x, y)

    plt.figure()
    plt.contourf(x, y, ux)
    plt.contour(x, y, ux, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.contourf(x, y, uy)
    plt.contour(x, y, uy, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.quiver(x, y, ux, uy)
    plt.ylim([-1.1,1.1])
    plt.show()

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
        solver_tol = 1e-6
    )
    return es, params

# def points():
#     points_grid([0, L, 20], [-c, c, 20], pts_filename)
#
def test_beam_bend():
    es, params = create_problem()
    problem = Controller(2, es, params)
    check_error(problem, 'traction', 'displacement', disp_bc, 4e-2)

    # # points()

    # interior_run(input_filename, pts_filename)
    # disp_intfilename = test_data_dir + 'beam_bend.disp_out_interior'
    # check_field(disp_intfilename, disp_bc, False, 6)

if __name__ == "__main__":
    plotter()
