import tbempy.TwoD
from elastic.compute import *
from elastic.input_builder import *

elements1 = [
    Element(
        pts = [[0, 0], [1, 0]],
        bc = [[0, 0], [0, 5]],
        bc_type = 'traction',
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
        pts = [[1, 0], [2, 0]],
        bc = [[1, 0], [1, 0]],
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

def test_form_linear_system():
    input = build_input(tbempy.TwoD, elements1, dict())
    systems = form_linear_systems(tbempy.TwoD, input)
    assert(len(systems[0]['lhs']) == 2)
    assert(len(systems[1]['lhs']) == 2)
    assert(systems[0]['rhs'].storage[0].storage[0] != 0.0)

def test_setup_integral_equation():
    input = build_input(tbempy.TwoD, elements1, dict())
    setup_bie = setup_integral_equation(tbempy.TwoD, input, input.bies[0])
    # There should be n_terms operators + one for the mass operator
    assert(len(setup_bie) == len(input.bies[0]['terms']) + 1)

def test_fields_from_bcs():
    meshes, bcs = meshes_bcs_from_elements(tbempy.TwoD, elements1)
    fields = fields_from_bcs(bcs)
    for k in bcs:
        assert((k, k) in fields)
        assert(len(fields[(k, k)].storage) == 2)
        assert(fields[(k, k)].storage[0].storage.shape[0] == bcs[k].size / 2)

#Regression test
def test_form_linear_system_no_traction_elements():
    input = build_input(tbempy.TwoD, elements2, dict())
    systems = form_linear_systems(tbempy.TwoD, input)
    # This segfaulted, so no assertions, just check if it runs

