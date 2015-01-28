import os
import h5py
import numpy as np
import subprocess
from mako.template import Template
import matplotlib.pyplot as plt

class Element(object):
    def __init__(self, pts, bc_type, bc, refine):
        self.pts = np.array(pts).astype(np.float32).tolist()
        self.bc_type = str(bc_type)
        self.bc = np.array(bc).astype(np.float32).tolist()
        self.refine = int(refine)

def line(end_pts, refine, bc_type, fnc):
    n = 2 ** refine
    x_vals = np.linspace(end_pts[0][0], end_pts[1][0], n + 1)
    y_vals = np.linspace(end_pts[0][1], end_pts[1][1], n + 1)
    ux, uy = fnc(x_vals, y_vals)

    es = []
    for i in range(n):
        es.append(Element(
            [[x_vals[i], y_vals[i]], [x_vals[i + 1], y_vals[i + 1]]],
            "displacement",
            [[ux[i], uy[i]], [ux[i + 1], uy[i + 1]]],
            0
        ))

    return es

def circle(center, r, refine, bc_type, fnc, reverse):
    n = 2 ** refine

    end_pt = 2 * np.pi
    if reverse:
        end_pt = -end_pt
    t = np.linspace(0.0, end_pt, n + 1)
    x = r * np.cos(t) + center[0]
    y = r * np.sin(t) + center[1]
    ux, uy = fnc(x, y)

    es = []
    for i in range(n):
        es.append(Element(
            [[x[i], y[i]], [x[i + 1], y[i + 1]]],
            bc_type,
            [[ux[i], uy[i]], [ux[i + 1], uy[i + 1]]],
            0
        ))

    return es


def exec_template(filename, **params):
    file_template = """
    {
        "shear_modulus": ${G},
        "poisson_ratio": ${mu},
        "elements": [
        % for e in es:
            {
                "pts": ${e.pts},
                "bc_type": "${e.bc_type}",
                "bc": ${e.bc},
                "refine": ${e.refine}
            }
            % if loop.index != len(es) - 1:
            ,
            % endif
        % endfor
        ]
    }
    """

    text = Template(file_template).render(**params)
    with open(filename, 'w') as file:
        file.write(text)

def run_file(filename, suppress_out):
    popen_params = dict(shell = True)
    if suppress_out:
        popen_params['stdout'] = subprocess.PIPE
    process = subprocess.Popen('./run ' + filename, **popen_params)
    process.wait()

def get_points(f):
    if f['locations'].shape[1] == 2:
        x = f['locations'][:, 0]
        y = f['locations'][:, 1]
        return x, y
    elif f['locations'].shape[1] == 4:
        x_index = [0, 2]
        y_index = [1, 3]
        vertices = np.array([
            f['locations'][:, x_index].flatten(),
            f['locations'][:, y_index].flatten()
        ]).T
        return vertices[:, 0], vertices[:, 1]

def check_field(filename, solution, plot_diff, digits):
    f = h5py.File(filename)

    x, y = get_points(f)
    datax = f['values0'][:, 0]
    datay = f['values1'][:, 0]
    exactx, exacty = solution(x, y)
    #TODO: It would be better to use norms here
    diffx = np.abs(exactx - datax)
    diffy = np.abs(exacty - datay)
    if plot_diff:
        plt.figure()
        plt.quiver(x, y, datax, datay)
        plt.figure()
        plt.quiver(x, y, exactx, exacty)
        plt.figure()
        plt.quiver(x, y, exactx - datax, exacty - datay)
        plt.show()
    np.testing.assert_almost_equal(diffx, np.zeros_like(diffx), digits)
    np.testing.assert_almost_equal(diffy, np.zeros_like(diffy), digits)
