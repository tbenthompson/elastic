import h5py
import numpy as np
import matplotlib.pyplot as plt
from input_builder import Element, line, exec_template, run_file, check_field
import subprocess

def G_from_E_mu(E, mu):
    return E / (2 * (1 + mu))

L = 1.0
c = 1.0
I = (2.0 / 3.0) * c ** 3

P = -40e6
E = 80e9
mu = 0.25
G = G_from_E_mu(E, mu)

# transformation from plane stress to plane strain
mu_fic = mu / (1 - mu)
E_fic = E / (1 - mu ** 2)
G_fic = G_from_E_mu(E_fic, mu_fic)

input_filename = 'test_data/beam_bend.in'

def disp_bc(x, y):
    ux = (-P * x ** 2 * y) / (2 * E_fic * I) \
          - (mu_fic * P * y ** 3) / (6 * E_fic * I) \
          + (P * y ** 3) / (6 * I * G_fic) \
          + (P * L ** 2 * y) / (2 * E_fic * I) \
          - (P * c ** 2 * y) / (2 * I * G_fic)

    uy = (mu_fic * P * x * y ** 2) / (2 * E_fic * I) \
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

def create_file():
    refine = 7

    es = []
    es.extend(line([[0, c], [0, -c]], refine, "displacement", disp_bc))
    es.append(Element([[0, -c], [L, -c]], "traction", [[0, 0], [0, 0]], refine))
    es.extend(line([[L, -c], [L, c]], refine, "displacement", disp_bc))
    es.append(Element([[L, c], [0, c]], "traction", [[0, 0], [0, 0]], refine))

    exec_template(input_filename, es = es, G = G, mu = mu)

def test_beam_bend():
    create_file()
    run_file(input_filename, stdout_dest = subprocess.PIPE)
    disp_filename = 'test_data/beam_bend.disp_out'
    check_field(disp_filename, disp_bc, False, 6)
    disp_intfilename = 'test_data/beam_bend.disp_outint'
    check_field(disp_intfilename, disp_bc, False, 6)

if __name__ == "__main__":
    plotter()
