import bie_spec
from collections import namedtuple
import numpy as np

Element = namedtuple('Element',
    ['pts', 'bc', 'bc_type', 'n_refines']
)
RefinedElement = namedtuple('RefinedElement',
    ['pts_mesh', 'bc_mesh', 'original_element']
)
#TODO: Remove all_mesh
Input = namedtuple('Input',
    [
        'elements', 'params', 'meshes', 'bcs', 'kernels',
        'quad_strategy', 'bies', 'all_mesh'
    ]
)

'''
Transforms input elements and parameters into meshes, boundary conditions,
kernels, quadrature formulae and boundary integral equations suitable for
input into tbem
'''
def build_input(tbem, elements, input_params):
    params = add_default_parameters(input_params)
    meshes, bcs = meshes_bcs_from_elements(tbem, elements)
    kernels = get_elastic_kernels(tbem, params)
    quad_strategy = get_quad_strategy(tbem, params)
    bies = bie_spec.get_all_BIEs(params)
    all_mesh = tbem.Mesh.create_union(meshes.values())
    # The gravity mesh should only include the outer boundaries of the
    # domain considered. Otherwise, the gravitational body force will be
    # effectively double counted.
    gravity_mesh = tbem.Mesh.create_union([meshes['traction'], meshes['displacement']])
    if params['gravity']:
        meshes['gravity'] = gravity_mesh
        bcs['gravity'] = np.ones((meshes['gravity'].n_facets(), tbem.dim, tbem.dim))
    return Input(
        elements, params, meshes, bcs, kernels, quad_strategy, bies, all_mesh
    )

def add_default_parameters(input_params):
    params = bie_spec.default_params()
    for k, v in input_params.iteritems():
        params[k] = v
    return params

'''
Builds mesh and boundary condition to tbem from the input elements
'''
def meshes_bcs_from_elements(tbem, elements):
    refined = refine_elements(tbem, elements)
    separated_by_bc = separate_elements_by_bc(refined)
    meshes, bcs = union_meshes_bcs(tbem, separated_by_bc)
    return add_empty_meshes(tbem, meshes, bcs)

'''
Refine the vertices and boundary conditions as specified in the input
'''
def refine_elements(tbem, elements):
    out_elements = []
    for e in elements:
        refined_pts = tbem.Mesh(np.array([e.pts])).refine_repeatedly(e.n_refines)
        refined_bc = tbem.Mesh(np.array([e.bc])).refine_repeatedly(e.n_refines)
        out_elements.append(RefinedElement(
            pts_mesh = refined_pts,
            bc_mesh = refined_bc,
            original_element = e
        ))
    return out_elements

'''
Separates the input elements into groups by the boundary condition type.
'''
def separate_elements_by_bc(elements):
    separated_by_bc = dict()
    for e in elements:
        bc_type = e.original_element.bc_type
        if bc_type in separated_by_bc:
            separated_by_bc[bc_type].append(e)
        else:
            separated_by_bc[bc_type] = [e]
    return separated_by_bc

'''
For each mesh type, this function combines each portion of the mesh
into one mesh object.
'''
def union_meshes_bcs(tbem, separated_by_bc):
    meshes = dict()
    bcs = dict()
    for k in separated_by_bc:
        pts_meshes = [e.pts_mesh for e in separated_by_bc[k]]
        meshes[k] = tbem.Mesh.create_union(pts_meshes)
        bcs[k] = np.vstack([e.bc_mesh.facets for e in separated_by_bc[k]])
    return meshes, bcs

'''
Add an empty mesh for each mesh type that is not present in the input.
'''
def add_empty_meshes(tbem, meshes, bcs):
    for k in bie_spec.mesh_types():
        if k not in meshes:
            meshes[k] = tbem.Mesh(np.empty((0, tbem.dim, tbem.dim)))
    for k in bie_spec.bc_types():
        if k not in bcs:
            bcs[k] = np.empty((0, tbem.dim, tbem.dim))
    return meshes, bcs

'''
Returns the set of four elastic kernels (also called Green's functions or
fundamental solutions) for the Somigliana identity and the hypersingular
integral equation
'''
def get_elastic_kernels(tbem, params):
    mu = params['shear_modulus']
    pr = params['poisson_ratio']
    kernels = dict(
        displacement = tbem.ElasticDisplacement(mu, pr),
        traction = tbem.ElasticTraction(mu, pr),
        adjoint_traction = tbem.ElasticAdjointTraction(mu, pr),
        hypersingular = tbem.ElasticHypersingular(mu, pr),
    )
    if params['gravity']:
        g = params['gravity_vector']
        kernels.update(dict(
            gravity_displacement = tbem.GravityDisplacement(mu, pr, g),
            gravity_traction = tbem.GravityTraction(mu, pr, g)
        ))
    return kernels

def get_quad_strategy(tbem, params):
    return tbem.QuadStrategy(
        params['obs_order'],
        params['singular_steps'],
        params['far_threshold'],
        1e-4
    )
