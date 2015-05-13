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
    assert(len(systems[0]['lhs']) == 4)
    assert(len(systems[1]['lhs']) == 4)
    assert(systems[0]['rhs'][0] != 0.0)

def test_setup_integral_equation():
    input = build_input(tbempy.TwoD, elements1, dict())
    setup_bie = setup_integral_equation(tbempy.TwoD, input, input.bies[0])
    # There should be n_terms operators + one for the mass operator
    assert(len(setup_bie) == len(input.bies[0]['terms']) + 1)

def test_fields_from_bcs():
    bcs = dict(abc = np.array([[[0, 1], [2, 3]], [[4, 5], [6, 7]]]))
    fields = fields_from_bcs(bcs)
    correct = np.array([[0, 2, 4, 6], [1, 3, 5, 7]])
    test = np.array(fields[('abc', 'abc')])
    assert(np.all(test == correct))

def test_real_fields_from_bcs():
    meshes, bcs = meshes_bcs_from_elements(tbempy.TwoD, elements1)
    fields = fields_from_bcs(bcs)
    for k in bcs:
        assert((k, k) in fields)
        n_components = len(fields[(k, k)])
        assert(n_components == 2)
        n_vals = fields[(k, k)][0].shape[0]
        assert(n_vals == bcs[k].size / n_components)

#Regression test
def test_form_linear_system_no_traction_elements():
    input = build_input(tbempy.TwoD, elements2, dict())
    systems = form_linear_systems(tbempy.TwoD, input)
    # This segfaulted, so no assertions, just check if it runs

