import h5py
import numpy as np
import matplotlib.pyplot as plt
from input_builder import Element, displacement_edge, exec_template, run_file

def G_from_E_mu(E, mu):
    return E / (2 * (1 + mu))

L = 1.0
c = 1.0
I = (2.0 / 3.0) * c ** 3

P = -10e6
E = 80e9
mu = 0.25
G = G_from_E_mu(E, mu)

# transformation from plane stress to plane strain
mu_fic = mu / (1 - mu)
E_fic = E / (1 - mu ** 2)
G_fic = G_from_E_mu(E_fic, mu_fic)

input_filename = 'test_data/beam_bend.in'

def solution(x, y):
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

def plotter():
    nx = 15
    ny = 15
    x_vals = np.linspace(0.0, L, nx)
    y_vals = np.linspace(-c, c, ny)
    x, y = np.meshgrid(x_vals, y_vals)

    ux, uy = solution(x, y)

    plt.figure()
    plt.contourf(x, y, ux)
    plt.contour(x, y, ux, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.contourf(x, y, uy)
    plt.contour(x, y, uy, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.quiver(x, y, ux, uy)
    plt.show()

def create_file():
    refine = 7
    # top_edge =
    # left_edge = Element([[0, -c], [0, c]], "traction", [[0, P], [0, P]], refine)

    es = []
    es.extend(displacement_edge([[0, c], [0, -c]], refine, solution))
    es.extend(displacement_edge([[0, -c], [L, -c]], refine, solution))
    # es.append(Element([[0, -c], [L, -c]], "traction", [[0, 0], [0, 0]], refine))
    es.extend(displacement_edge([[L, -c], [L, c]], refine, solution))
    es.extend(displacement_edge([[L, c], [0, c]], refine, solution))
    # es.append(Element([[L, c], [0, c]], "traction", [[0, 0], [0, 0]], refine))

    # fictitious_mu = mu / (1 - mu)
    # fictitious_E = E / (1 - mu ** 2)
    # fictitious_G = G_from_E_mu(fictitious_E, fictitious_mu)
    # exec_template(input_filename, es = es, G = fictitious_G, mu = fictitious_mu)
    exec_template(input_filename, es = es, G = G, mu = mu)
    print("Input file created")

def test_displacements():
    filename = 'test_data/beam_bend.disp_outint'
    f = h5py.File(filename)
    x = f['locations'][:, 0]
    y = f['locations'][:, 1]
    datax = f['values0'][:, 0]
    datay = f['values1'][:, 0]
    exactx, exacty = solution(x, y)
    errorx = np.abs((exactx - datax) / exactx)
    errory = np.abs((exacty - datay) / exacty)
    np.testing.assert_almost_equal(errorx, np.zeros_like(errorx), 2)
    np.testing.assert_almost_equal(errory, np.zeros_like(errory), 2)
    print("Tests passed!")


if __name__ == "__main__":
    # plotter()
    create_file()
    run_file(input_filename)
    test_displacements()

