import subprocess
import os
from input_builder import Element, bem_template, run, check_field
import numpy as np

def build_slip_bc(a, stress_drop, G, mu):
    def slip_bc(x, y):
        # Classic Griffith crack solution!
        slip = 2 * (1 - mu) * (stress_drop / G) * np.sqrt(a ** 2 - x ** 2)
        return slip, np.zeros_like(x)
    return slip_bc

def create_file(a, stress_drop, G, mu, input_filename):
    refine = 7
    es = []
    es.append(Element([[-a, 0], [-a/1.1, 0]], "crack_traction",
        [[stress_drop, 0], [stress_drop, 0]], refine))
    es.append(Element([[-a/1.1, 0], [a/1.1, 0]], "crack_traction",
        [[stress_drop, 0], [stress_drop, 0]], refine))
    es.append(Element([[a/1.1, 0], [a, 0]], "crack_traction",
        [[stress_drop, 0], [stress_drop, 0]], refine))
    bem_template(input_filename, es = es, G = G, mu = mu)

def test_crack():
    input_filename = 'test_data/crack.in'
    a = 10000
    G = 40e9
    mu = 0.29
    stress_drop = 5.5e6

    in_root, file_ext = os.path.splitext(input_filename)
    slip_filename = in_root + '.slip_out'

    slip_bc = build_slip_bc(a, stress_drop, G, mu)
    create_file(a, stress_drop, G, mu, input_filename)
    run(input_filename, stdout_dest = subprocess.PIPE)
    check_field(slip_filename, slip_bc, False, 5)
