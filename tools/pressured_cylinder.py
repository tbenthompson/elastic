import numpy as np
import matplotlib.pyplot as plt
from input_builder import circle, exec_template, run_file, test_displacements

a = 0.3
b = 1.2

p_a = 10e6
p_b = -50e6

E = 80e9
mu = 0.25
G = E / (2 * (1 + mu))

input_filename = 'test_data/pressured_cylinder.in'

def solution(x, y):
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    ur = ((1 + mu) * a ** 2 * b ** 2) / (E * (b ** 2 - a ** 2)) *\
        (((p_a - p_b) / r) +
         ((1 - 2 * mu) * ((p_a * a ** 2 - p_b * b ** 2) / (a ** 2 * b ** 2)) * r))
    ux = ur * np.cos(theta)
    uy = ur * np.sin(theta)

    return ux, uy

def trac_bc(pressure):
    def bc(x, y):
        theta = np.arctan2(y, x)
        p_x = pressure * np.cos(theta)
        p_y = pressure * np.sin(theta)
        return p_x, p_y
    return bc

def plotter():
    nt = 20
    nr = 5
    t_vals = np.linspace(0.0, 2 * np.pi, nt)
    r_vals = np.linspace(a, b, nr)
    r, t = np.meshgrid(r_vals, t_vals)
    x = r * np.cos(t)
    y = r * np.sin(t)

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

def create_file(bc_types):
    refine = 8

    es = []
    if bc_types["inner"] == "traction":
        es.extend(circle([0, 0], a, refine, "traction", trac_bc(-p_a), True))
    elif bc_types["inner"] == "displacement":
        es.extend(circle([0, 0], a, refine, "displacement", solution, True))

    if bc_types["outer"] == "traction":
        es.extend(circle([0, 0], b, refine, "traction", trac_bc(p_b), False))
    elif bc_types["outer"] == "displacement":
        es.extend(circle([0, 0], b, refine, "displacement", solution, False))

    exec_template(input_filename, es = es, G = G, mu = mu)
    print("Input file created")

if __name__ == "__main__":
    out_filename = 'test_data/pressured_cylinder.disp_outint'
    # plotter()

    create_file(dict(inner = "traction", outer = "traction"))
    run_file(input_filename, True)
    test_displacements(out_filename, solution, False)

    create_file(dict(inner = "displacement", outer = "displacement"))
    run_file(input_filename, True)
    test_displacements(out_filename, solution, False)

    # create_file(dict(inner = "traction", outer = "displacement"))
    # run_file(input_filename, True)
    # test_displacements(out_filename, solution, False)

    # create_file(dict(inner = "displacement", outer = "traction"))
    # run_file(input_filename, True)
    # test_displacements(out_filename, solution, False)
