import tbempy.TwoD
import numpy as np
from elastic.bie_spec import field_types, get_elastic_kernels, get_BIEs
from elastic.constraints import form_displacement_constraints, gather_bc_constraints,\
    Term, transform_element_local_term, any_displacement_constraints
from elastic.compute import IntegralDispatcher, Op
from elastic.system import split_into_components, scale
from elastic.dof_handling import DOFMap
from elastic.element_types import displacement, traction, slip
from elastic.meshing import build_meshes, line, postprocess_meshes

es = [
    displacement([[0, 0], [1, 0]], [[1, 2], [3, 4]]),
    traction([[1, 1], [2, 2]], [[5, 6], [7, 8]])
]
dof_map_internal = dict()
dof_map_internal[('continuous', 'displacement')] = [7, 13]
dof_map_internal[('continuous', 'traction')] = [20, 22]
dof_map = DOFMap(2, 24, dof_map_internal)
meshes = dict(
    continuous = tbempy.TwoD.line_mesh([0, 0], [10, 0]).refine_repeatedly(4),
)
kernels = get_elastic_kernels(
    tbempy.TwoD, dict(
        shear_modulus = 30e9, poisson_ratio = 0.25,
        gravity_vector = [0.0, -9.8 * 2700]
    )
)

def test_displacement_constraints_continuity():
    meshes['discontinuous'] = tbempy.TwoD.Mesh(np.empty((0, 2, 2)))
    result = form_displacement_constraints(tbempy.TwoD, dof_map, meshes)
    assert(len(result) == 2 * (meshes['continuous'].n_facets() - 1))

def test_displacement_constraints_fault():
    meshes['discontinuous'] = tbempy.TwoD.line_mesh([5, 0], [5, 5])
    result = form_displacement_constraints(tbempy.TwoD, dof_map, meshes)
    assert(len(result) == 2 * (meshes['continuous'].n_facets() - 2))

def test_build_dof_map():
    meshes['discontinuous'] = tbempy.TwoD.line_mesh([5, 0], [5, 5])
    dof_map = DOFMap.build(2, field_types, meshes)
    disp_dofs = dof_map.map[('continuous', 'displacement')]
    assert(disp_dofs[1] - disp_dofs[0] == meshes['continuous'].n_dofs())

def test_build_dof_map_past_end():
    dof_map = DOFMap.build(2, field_types, meshes)
    disp_dofs = dof_map.map[('continuous', 'displacement')]
    assert(disp_dofs[2] - disp_dofs[1] == meshes['continuous'].n_dofs())

def test_expand():
    dof_map = DOFMap(2, 4, {
        ('A', '1'): [0, 1, 2],
        ('B', '2'): [2, 3, 4]
    })
    result = dof_map.expand([1,2,3,4])
    np.testing.assert_equal(result[('A', '1')], [[1], [2]])
    np.testing.assert_equal(result[('B', '2')], [[3], [4]])

def test_concatenate():
    dof_map = DOFMap(2, 4, {
        ('A', '1'): [0, 1, 2],
        ('B', '2'): [2, 3, 4]
    })
    input = {
        ('A', '1'): [[0], [1]],
        ('B', '2'): [[2], [3]]
    }
    np.testing.assert_equal(dof_map.concatenate(input), [0, 1, 2, 3])

def test_build_meshes():
    result, element_lists = build_meshes(tbempy.TwoD, es)
    assert(result['continuous'].n_facets() == 2)
    assert(result['discontinuous'].n_facets() == 0)

def test_gather_bc_constraints():
    result = gather_bc_constraints(tbempy.TwoD, dof_map, es)
    correct = [
        (7, 1.0, 1.0), (13, 1.0, 2.0), (8, 1.0, 3.0), (14, 1.0, 4.0),
        (22, 1.0, 5.0), (24, 1.0, 6.0), (23, 1.0, 7.0), (25, 1.0, 8.0)
    ]
    for r, c in zip(result, correct):
        assert(r.terms[0].dof == c[0])
        assert(r.terms[0].weight == c[1])
        assert(r.rhs == c[2])

class FakeEvaluator(object):
    def __init__(self):
        self.record = []
    def mass(self, obs_mesh):
        self.record.append(obs_mesh)
    def boundary(self, obs_mesh, src_mesh, kernel):
        self.record.append((obs_mesh, src_mesh, kernel))
        return self.record[-1]
    def interior(self, pts, normals, src_mesh, kernel):
        self.record.append((pts, normals, src_mesh, kernel))
        return self.record[-1]

def test_integral_dispatcher():
    evaluator = FakeEvaluator()
    dispatcher = IntegralDispatcher(meshes, kernels, evaluator)
    dispatcher.compute_boundary(dict(
        obs_mesh = 'continuous', src_mesh = 'discontinuous', kernel = 'hypersingular'
    ))
    dispatcher.compute_interior(
        dict(src_mesh = 'discontinuous', kernel = 'traction'),
        'pts', 'normals'
    )
    dispatcher.compute_mass(dict(src_mesh = 'continuous'))

def test_line():
    mesh = line(
        [[0, 0], [1, 0]], 3,
        lambda pts: displacement(pts, [[0,0],[0,0]])
    )
    assert(mesh[0]['type'] == 'continuous')
    np.testing.assert_almost_equal(mesh[0]['pts'], [[0, 0], [0.125, 0]])

def test_get_bies():
    bies = get_BIEs(dict(
        gravity = False, shear_modulus = 1, poisson_ratio = 0.25,
        length_scale = 1
    ))

def scale_tester(inverse):
    unknowns = dict()
    unknowns[('A', '1')] = [4, 5, 6]
    unknowns[('B', '2')] = [8, 10, 12]
    field_types = dict()
    field_types[('A', '1')] = lambda p: 2.0
    field_types[('B', '2')] = lambda p: 1.0
    if inverse:
        field_types[('A', '1')] = lambda p: 1.0
        field_types[('B', '2')] = lambda p: 2.0
    scale(unknowns, field_types, None, inverse)
    np.testing.assert_equal(unknowns.values()[0], unknowns.values()[1])

def test_scale():
    scale_tester(False)

def test_scale_inverse():
    scale_tester(True)

def test_op_apply():
    class Fake(object):
        def apply(self, stuff):
            return 1.0
    thing = Op(Fake(), dict(multiplier = 3.0), False)
    assert(thing.apply([np.array([0])]) == 3.0)

def test_split_into_components():
    result = split_into_components(2, np.array([1,2,3,4]))
    np.testing.assert_equal(result, [[1,2],[3,4]])

def test_transform_element_local_term_normal():
    t = Term('slip', 'normal', 1.0)
    result = transform_element_local_term(2, t, [[0, 0], [1, 0]])
    assert(result[0].component == 0)
    assert(result[0].weight == 0)
    assert(result[1].component == 1)
    assert(result[1].weight == 1)

def test_transform_element_local_term_multiplier():
    t = Term('slip', 'normal', -2.0)
    result = transform_element_local_term(2, t, [[0, 0], [1, 0]])
    # assert(result[1].component ==

def test_transform_element_local_term_tangential():
    t = Term('slip', 'tangential0', 1.0)
    result = transform_element_local_term(2, t, [[0, 0], [1, 0]])

def test_build_meshes_cut_at_fault():
    meshes = postprocess_meshes(tbempy.TwoD, build_meshes(tbempy.TwoD, [
        dict(type = 'continuous', pts = [[-1, 0], [1, 0]]),
        dict(type = 'discontinuous', pts = [[0, -1], [0, 1]])
    ])[0])
    assert(meshes['continuous'].n_facets() == 2)
    assert(meshes['discontinuous'].n_facets() == 2)

def test_any_displacement_constraints():
    es_no_disp = [
        traction([[1, 1], [2, 2]], [[5, 6], [7, 8]])
    ]
    assert(any_displacement_constraints(es))
    assert(not any_displacement_constraints(es_no_disp))
