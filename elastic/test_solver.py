import tbempy.TwoD
from elastic.solver import *
from elastic.compute import form_linear_systems
from elastic.input_builder import *

elements1 = [
    Element(
        pts = [[0, 0], [1, 0]],
        bc = [[5, 0], [0, 0]],
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

def test_dof_map():
    input = build_input(tbempy.TwoD, elements1, dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    correct_starts = [0, 2, 4, 6]
    assert(np.all(dof_map.start_positions == correct_starts))
    assert(dof_map.n_components == 4)
    assert(dof_map.n_dofs == 8)

def test_constraint_matrix():
    input = build_input(tbempy.TwoD, elements1, dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    constraint_matrix = build_constraint_matrix(
        tbempy.TwoD, dof_map, input.bies, input.meshes
    )
    # Nothing from constraint_matrix is exposed to python so no tests are
    # performed except to check that the creation succeeds
    assert(constraint_matrix is not None)

def test_concatenate_condense():
    input = build_input(tbempy.TwoD, elements1, dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    constraint_matrix = build_constraint_matrix(
        tbempy.TwoD, dof_map, input.bies, input.meshes
    )
    systems = form_linear_systems(tbempy.TwoD, input)
    result = concatenate_condense(
        dof_map, constraint_matrix, [s['rhs'] for s in systems]
    )
