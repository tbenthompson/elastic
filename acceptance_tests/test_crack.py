import subprocess
import os
import numpy as np
from tools.input_builder import *

def build_slip_bc(a, stress_drop, G, mu):
    def slip_bc(x, y):
        # Classic Griffith crack solution! Normally, the solution is given
        # as a displacement for each crack face. I multiply by 2 to get the
        # slip across the slip.
        slip = 2 * (1 - mu) * (stress_drop / G) * np.sqrt(a ** 2 - x ** 2)
        return slip, np.zeros_like(x)
    return slip_bc

def create_file(a, stress_drop, G, mu, input_filename, bc_type):
    refine = 7
    es = []
    es.append(Element([[-a, 0], [-a/1.1, 0]], bc_type,
        [[stress_drop, 0], [stress_drop, 0]], refine))
    es.append(Element([[-a/1.1, 0], [a/1.1, 0]], bc_type,
        [[stress_drop, 0], [stress_drop, 0]], refine))
    es.append(Element([[a/1.1, 0], [a, 0]], bc_type,
        [[stress_drop, 0], [stress_drop, 0]], refine))
    bem_template(input_filename, es, shear_modulus = G, mu = mu)

def griffith_soln(bc_type):
    input_filename = test_data_dir + bc_type + '_griffith.in'
    a = 0.5
    G = 40e9
    mu = 0.29
    stress_drop = 5.5e6

    in_root, file_ext = os.path.splitext(input_filename)
    slip_filename = in_root + '.slip_out'

    slip_bc = build_slip_bc(a, stress_drop, G, mu)
    create_file(a, stress_drop, G, mu, input_filename, bc_type)
    run(input_filename, stdout_dest = subprocess.PIPE)
    check_field(slip_filename, slip_bc, False, 5)

def test_crack():
    griffith_soln("crack_traction")

# def test_free_slip():
#     griffith_soln("free_slip")
