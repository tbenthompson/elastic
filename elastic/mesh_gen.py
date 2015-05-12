import numpy as np
from elastic.input_builder import Element

#TODO: Lots of overlap with tbempy mesh_gen
def line(end_pts, refine, bc_type, fnc):
    n = 2 ** refine
    x_vals = np.linspace(end_pts[0][0], end_pts[1][0], n + 1)
    y_vals = np.linspace(end_pts[0][1], end_pts[1][1], n + 1)
    ux, uy = fnc([x_vals, y_vals])

    es = []
    for i in range(n):
        #TODO: Why hard coded as "displacement"?
        es.append(Element(
            [[x_vals[i], y_vals[i]], [x_vals[i + 1], y_vals[i + 1]]],
            [[ux[i], uy[i]], [ux[i + 1], uy[i + 1]]],
            bc_type,
            0
        ))

    return es

def circle(center, r, refine, bc_type, fnc, reverse):
    # If refine is 0, then by adding two, the approximation is still
    # non-degenerate
    n = 2 ** (refine + 2)


    # To reverse the orientation of the circle, simply go from
    # 0 --> -2pi rather than 0 --> 2pi
    end_pt = 2 * np.pi
    if reverse:
        end_pt = -end_pt

    t = np.linspace(0.0, end_pt, n + 1)
    x = r * np.cos(t) + center[0]
    y = r * np.sin(t) + center[1]

    ux, uy = fnc([x, y])

    es = []
    for i in range(n):
        vs = [[x[i], y[i]], [x[i + 1], y[i + 1]]]
        es.append(Element(
            vs,
            [[ux[i], uy[i]], [ux[i + 1], uy[i + 1]]],
            bc_type,
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
        bcs = [fnc(v) for v in vs]
        es.append(Element(vs, bcs, bc_type, 0))
    return es
