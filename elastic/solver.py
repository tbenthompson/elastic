from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD
import bie_spec
from collections import namedtuple
import numpy as np

Element = namedtuple('Element',
    ['pts', 'bc', 'bc_type', 'n_refines']
)
RefinedElement = namedtuple('RefinedElement',
    ['pts_mesh', 'bc_mesh', 'original_element']
)

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def bc_types():
    return [
        'traction',
        'displacement',
        'slip',
        'crack_traction',
        'free_slip_traction'
    ]

def mesh_types():
    return bc_types()

def meshes_from_elements(tbem, elements):
    refined = refine_elements(elements)
    separated_by_bc = separate_elements_by_bc(refined)
    meshes = dict()
    bcs = dict()
    for k in separated_by_bc:
        pts_meshes = [e.pts_mesh for e in separated_by_bc[k]]
        bc_meshes = [e.bc_mesh for e in separated_by_bc[k]]
        meshes[k] = tbem.Mesh.create_union(pts_meshes)
    for k in mesh_types():
        if k not in meshes:
            meshes[k] = tbem.Mesh(np.empty((0, 2, 2)))
    return meshes

def refine_elements(elements):
    out_elements = []
    for e in out_elements:
        refined_pts = tbem.Mesh(np.array([e.pts])).refine_repeatedly(e.n_refines)
        refined_bc = tbem.Mesh(np.array([e.bc])).refine_repeatedly(e.n_refines)
        out_elements.append(RefinedElement(
            pts_mesh = refined_pts,
            bc_mesh = refined_bc,
            original_element = e
        ))
    return out_elements

def separate_elements_by_bc(elements):
    separated_by_bc = dict()
    for e in elements:
        bc_type = e.original_element.bc_type
        if bc_type in separated_by_bc:
            separated_by_bc[bc_type].append(e)
        else:
            separated_by_bc[bc_type] = [e]
    return separated_by_bc

def default_params():
    return dict(
        obs_order = 3,
        singular_steps = 8,
        far_threshold = 3.0,
        near_tol = 1e-4,
        solver_tol = 1e-6,
        poisson_ratio = 0.25,
        shear_modulus = 30e9
    )

def solve(dim, elements, params):
    tbem = get_tbem(dim)
    meshes = meshes_from_elements(tbem, elements)
    for k,v in meshes.iteritems():
        print(k + ': ' + str(v.facets))
    bies = bie_spec.get_all_BIEs()
    dof_map = build_block_dof_map([2, 3])

if __name__ == "__main__":
    solve(2, [
        Element(
            pts = [[0, 0], [1, 0]],
            bc = [[0, 0], [0, 0]],
            bc_type = 'displacement',
            n_refines = 0
        )
    ], {})

