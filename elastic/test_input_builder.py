import tbempy.TwoD
from elastic.input_builder import *
from elastic.bie_spec import *

elements1 = [
    Element(
        pts = [[0, 0], [1, 0]],
        bc = [[0, 0], [0, 0]],
        bc_type = 'displacement',
        n_refines = 0
    ),
    Element(
        pts = [[1, 0], [2, 0]],
        bc = [[1, 0], [1, 0]],
        bc_type = 'displacement',
        n_refines = 0
    )
]

elements2 = [
    Element(
        pts = [[0, 0], [1, 0]],
        bc = [[0, 0], [0, 0]],
        bc_type = 'displacement',
        n_refines = 1
    ),
    Element(
        pts = [[1, 0], [2, 0]],
        bc = [[1, 0], [2, 0]],
        bc_type = 'displacement',
        n_refines = 1
    )
]

params1 = dict(shear_modulus = 10e9, singular_steps = 10)

def test_meshes_bcs_from_elements():
    meshes, bcs = meshes_bcs_from_elements(tbempy.TwoD, elements1)
    assert(meshes['displacement'].facets.shape[0] == 2)
    assert(bcs['displacement'].shape[0] == 2)

def test_meshes_bcs_from_elements_all_meshes():
    meshes, bcs = meshes_bcs_from_elements(tbempy.TwoD, elements1)
    for m in mesh_types():
        assert(m in meshes)
    for m in bc_types():
        assert(m in bcs)

def test_meshes_bcs_from_elements_refine():
    meshes, bcs = meshes_bcs_from_elements(tbempy.TwoD, elements2)
    assert(meshes['displacement'].facets.shape[0] == 4)
    assert(bcs['displacement'].shape[0] == 4)
    assert(bcs['displacement'][2][1][0] == 1.5)

def test_elastic_kernels():
    p = dict(shear_modulus = 30e9, poisson_ratio = 0.25)
    kernels = get_elastic_kernels(tbempy.TwoD, p)
    assert(len(kernels.keys()) == 4)

def test_input_parameters():
    p = add_default_parameters(params1)
    assert(p['obs_order'] == 3)
    assert(p['singular_steps'] == 10)
    assert(p['shear_modulus'] == 10e9)

def test_full_input_build():
    input = build_input(tbempy.TwoD, elements2, params1)
    assert(len(input.kernels.keys()) == 4)
    assert(len(input.bies) == 2)
