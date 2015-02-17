import os
import h5py
import numpy as np
import subprocess
from mako.template import Template
import matplotlib.pyplot as plt

test_data_dir = 'test_data/auto_gen/'

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

####TODO: This code is duplicated with mesh_gen.cpp in the library
def sphere_midpt(pt1, pt2, center, r):
    midpt = (np.array(pt1) + np.array(pt2)) / 2
    return (r / np.linalg.norm(midpt - center)) * (midpt - center) + center

def sphere(center, r, refine, bc_type, fnc, reverse):
    verts = [
        [0.0, -r, 0.0], [r, 0.0, 0.0], [0.0, 0.0, r],
        [-r, 0.0, 0.0], [0.0, 0.0, -r], [0.0, r, 0.0]
    ]
    faces = [
        [1, 0, 2], [2, 0, 3], [3, 0, 4], [4, 0, 1],
        [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]
    ]

    ## TODO: This code refines and is duplicated in mesh.cpp
    ## TODO: THEME: Need some sort of lightweight python-cpp integration
    for i in range(refine):
        new_faces = []
        for f in faces:
            midpts = [
                sphere_midpt(verts[f[0]], verts[f[1]], center, r),
                sphere_midpt(verts[f[1]], verts[f[2]], center, r),
                sphere_midpt(verts[f[2]], verts[f[0]], center, r)
            ]
            midptidx = [len(verts) + i for i in range(3)]
            verts.extend(midpts)
            new_faces.append([f[0], midptidx[0], midptidx[2]])
            new_faces.append([f[1], midptidx[1], midptidx[0]])
            new_faces.append([f[2], midptidx[2], midptidx[1]])
            new_faces.append([midptidx[0], midptidx[1], midptidx[2]])
        faces = new_faces

    es = []
    for f in faces:
        vs = [verts[f[i]] for i in range(3)]
        if reverse:
            vs = [vs[1], vs[0], vs[2]]
        bcs = [fnc(*v) for v in vs]
        es.append(Element(vs, bc_type, bcs, 0))
    return es

def make_points(x_range, y_range):
    x_vals = np.linspace(*x_range)
    y_vals = np.linspace(*y_range)
    x_pts, y_pts = np.meshgrid(x_vals, y_vals)
    ps = zip(x_pts.flatten(), y_pts.flatten())
    return ps

def points_grid(x_range, y_range, out_filename):
    points_template(out_filename, make_points(x_range, y_range))

def points_template(out_filename, ps):
    file_template = """
    [
        % for p in ps:
        [${p[0]}, ${p[1]}]
        % if loop.index != len(ps) - 1:
        ,
        % endif
        % endfor
    ]
    """
    exec_template(file_template, out_filename, ps = ps)

def bem_template(filename, es, **params):
    file_template = """
    {
        % for key,value in params.iteritems():
            "${key}": ${value},
        % endfor
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
    exec_template(file_template, filename, es = es, params = params)

def exec_template(file_template, filename, **params):
    text = Template(file_template).render(**params)
    with open(filename, 'w') as file:
        file.write(text)


def get_config():
    config = dict()
    config['solver_2d'] = './solve2d'
    config['solver_3d'] = './solve3d'
    config['interior_2d'] = './interior'
    return config

def run(filename, dim = 2, stdout_dest = None):
    execute(get_config()['solver_' + str(dim) + 'd'] + ' ' + filename, stdout_dest)

def interior_run(bem_filename, pts_filename, stdout_dest = None):
    execute(get_config()['interior_2d'] + ' ' + bem_filename + ' ' + pts_filename, stdout_dest)

def execute(cmd, stdout_dest):
    popen_params = dict(shell = True)
    popen_params['stdout'] = stdout_dest
    process = subprocess.Popen(cmd, **popen_params)
    process.wait()

def get_points(f):
    if f['locations'].shape[1] == 2:
        x = f['locations'][:, 0]
        y = f['locations'][:, 1]
        return x, y
    elif f['locations'].shape[1] == 3:
        x = f['locations'][:, 0]
        y = f['locations'][:, 1]
        z = f['locations'][:, 2]
        return x, y, z
    elif f['locations'].shape[1] == 4:
        x_index = [0, 2]
        y_index = [1, 3]
        vertices = np.array([
            f['locations'][:, x_index].flatten(),
            f['locations'][:, y_index].flatten()
        ]).T
        return vertices[:, 0], vertices[:, 1]
    elif f['locations'].shape[1] == 9:
        x_index = [0, 3, 6]
        y_index = [1, 4, 7]
        z_index = [2, 5, 8]
        vertices = np.array([
            f['locations'][:, x_index].flatten(),
            f['locations'][:, y_index].flatten(),
            f['locations'][:, z_index].flatten()
        ]).T
        return vertices[:, 0], vertices[:, 1], vertices[:, 2]

def check_field(filename, solution, plot_diff, digits, point_limiter = None):
    if point_limiter is None:
        point_limiter = lambda x, y: np.ones_like(x, dtype = np.bool)

    f = h5py.File(filename)

    pts = get_points(f)
    data = []
    exact = []
    diff = []
    exact = solution(*pts)
    for d in range(len(pts)):
        data = f['values' + str(d)][:, 0]
        diff = np.abs(exact[d] - data)
        np.testing.assert_almost_equal(diff, np.zeros_like(diff), digits)
