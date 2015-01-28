import numpy as np
import matplotlib.pyplot as plt
import os
from input_builder import circle, exec_template, run_file, check_field

def build_disp_bc(a, b, p_a, p_b, E, mu):
    def disp_bc(x, y):
        r = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan2(y, x)
        ur = ((1 + mu) * a ** 2 * b ** 2) / (E * (b ** 2 - a ** 2)) *\
            (((p_a - p_b) / r) +
             ((1 - 2 * mu) * ((p_a * a ** 2 - p_b * b ** 2) / (a ** 2 * b ** 2)) * r))
        ux = ur * np.cos(theta)
        uy = ur * np.sin(theta)
        return ux, uy
    return disp_bc

def build_trac_bc(a, b, p_a, p_b, E, mu):
    def trac_bc(x, y):
        theta = np.arctan2(y, x)
        r = np.sqrt(x ** 2 + y ** 2)
        pressure = np.where(r < ((a + b) / 2), -p_a, p_b)
        p_x = pressure * np.cos(theta)
        p_y = pressure * np.sin(theta)
        return p_x, p_y
    return trac_bc

def plotter(disp_bc):
    nt = 20
    nr = 5
    t_vals = np.linspace(0.0, 2 * np.pi, nt)
    r_vals = np.linspace(a, b, nr)
    r, t = np.meshgrid(r_vals, t_vals)
    x = r * np.cos(t)
    y = r * np.sin(t)

    ux, uy = disp_bc(x, y)

    plt.figure()
    plt.contourf(x, y, ux)
    plt.contour(x, y, ux, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.contourf(x, y, uy)
    plt.contour(x, y, uy, linestyles = 'solid', colors = 'k', linewidths = 2)
    plt.figure()
    plt.quiver(x, y, ux, uy)
    plt.show()

def create_file(a, b, E, mu, input_filename, bc_types, bc_funcs):
    refine = 9

    es = []
    es.extend(circle([0, 0], a, refine, bc_types['inner'],
                     bc_funcs[bc_types['inner']], True))
    es.extend(circle([0, 0], b, refine, bc_types['outer'],
                     bc_funcs[bc_types['outer']], False))

    G = E / (2 * (1 + mu))
    exec_template(input_filename, es = es, G = G, mu = mu)

def delete_files(input_filepath):
    dir, filename = input_filepath.split('/')
    root = os.path.splitext(filename)[0]
    for f in os.listdir(dir):
        if root in f:
            os.remove(os.path.join(dir, f))

def pressured_cylinder(bc_types):
    a = 0.3
    b = 1.2
    p_a = 10e6
    p_b = -15e6
    E = 80e9
    mu = 0.25
    input_filename = 'test_data/pressured_cylinder.in'

    disp_bc = build_disp_bc(a, b, p_a, p_b, E, mu)
    trac_bc = build_trac_bc(a, b, p_a, p_b, E, mu)
    bc_funcs = dict()
    bc_funcs['traction'] = trac_bc
    bc_funcs['displacement'] = disp_bc

    in_root, file_ext = os.path.splitext(input_filename)
    displacement_filename = in_root + '.disp_out'
    traction_filename = in_root + '.trac_out'
    interior_disp_filename = in_root + '.disp_outint'

    delete_files(input_filename)
    create_file(a, b, E, mu, input_filename,
                bc_types, bc_funcs)
    run_file(input_filename, False)
    if 'displacement' in bc_types.values():
        check_field(traction_filename, trac_bc, False, -5)
    if 'traction' in bc_types.values():
        check_field(displacement_filename, disp_bc, False, 7)
    # check_field(interior_disp_filename, disp_bc, False, 7)

def test_trac_trac():
    pressured_cylinder(dict(inner = "traction", outer = "traction"))

def test_disp_trac():
    pressured_cylinder(dict(inner = "displacement", outer = "traction"))

def test_trac_disp():
    pressured_cylinder(dict(inner = "traction", outer = "displacement"))

def test_disp_disp():
    pressured_cylinder(dict(inner = "displacement", outer = "displacement"))

if __name__ == "__main__":
    test_trac_disp()
