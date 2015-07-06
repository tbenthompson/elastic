import tbempy.TwoD
from elastic.dof_handling import *
from elastic.interface import *
from elastic.compute import form_linear_systems
from elastic.input_builder import *

def get_element_list():
    return [
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
    input = build_input(tbempy.TwoD, get_element_list(), dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    assert(dof_map[('traction', 'displacement')][0] == 4)
    assert(dof_map[('traction', 'displacement')][1] == 6)
    assert(dof_map[('traction', 'displacement')][2] == 8)


#TODO: These tests need some kind of fixture-based design so that
# they stop replicating the same data
def test_constraint_matrix():
    input = build_input(tbempy.TwoD, get_element_list(), dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    constraint_matrix = build_constraint_matrix(
        tbempy.TwoD, dof_map, input
    )
    #TODO: Check the values in the constraint matrix
    assert(constraint_matrix is not None)

def test_concatenate_condense():
    input = build_input(tbempy.TwoD, get_element_list(), dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    constraint_matrix = build_constraint_matrix(
        tbempy.TwoD, dof_map, input
    )
    systems = form_linear_systems(tbempy.TwoD, input)
    result = concatenate_condense(
        tbempy.TwoD, dof_map, constraint_matrix, [s['rhs'] for s in systems]
    )

def test_distribute_expand():
    elements = get_element_list()
    elements[1] = elements[1]._replace(bc_type = 'traction')
    input = build_input(tbempy.TwoD, elements, dict())
    dof_map = build_dof_map(tbempy.TwoD, input.bies, input.meshes)
    constraint_matrix = build_constraint_matrix(
        tbempy.TwoD, dof_map, input
    )
    systems = form_linear_systems(tbempy.TwoD, input)
    a = [0, 1, 0, 1, 0, 1]
    fields = distribute_expand(tbempy.TwoD, dof_map, constraint_matrix, a)

    assert(('displacement', 'traction') in fields)
    assert(fields[('displacement', 'traction')][0].shape[0] == 0)
    assert(fields[('displacement', 'traction')][1].shape[0] == 0)
    assert(('traction', 'displacement') in fields)
    assert(np.all(fields[('traction', 'displacement')][0] == [0, 1, 1, 0]))
    assert(np.all(fields[('traction', 'displacement')][1] == [1, 0, 0, 1]))
