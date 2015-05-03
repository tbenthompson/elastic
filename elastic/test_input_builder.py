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
        pts = [[0, 0], [1, 0]],
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
        pts = [[0, 0], [1, 0]],
        bc = [[1, 0], [2, 0]],
        bc_type = 'displacement',
        n_refines = 1
    )
]

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
    kernels = get_elastic_kernels(tbempy.TwoD, 30e9, 0.25)
    assert(len(kernels.keys()) == 4)

def test_input_parameters():
    p = form_parameters(dict(shear_modulus = 10e9, singular_steps = 10))
    assert(p['obs_order'] == 3)
    assert(p['singular_steps'] == 10)
    assert(p['shear_modulus'] == 10e9)

