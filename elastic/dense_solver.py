import numpy as np
from dof_handling import distribute_expand, concatenate_condense, \
    scale_columns, scale_rows, operate_on_solution_fields
import iterative_solver

def dense_solver(tbem, input, dof_map, constraint_matrix, systems):
    np_op = dense_matrix(tbem, input, dof_map, constraint_matrix, systems)
    rhs = dense_rhs(tbem, input, dof_map, constraint_matrix, systems)

    mat_vec = iterative_solver.build_matrix_vector_product(
        tbem, input, dof_map, constraint_matrix, systems
    )
    rhs -= mat_vec(np.zeros_like(rhs))

    input.logger.linear_solve_start(rhs.shape[0])
    soln = np.linalg.solve(np_op, rhs)
    input.logger.linear_solve_end()
    unknowns = distribute_expand(tbem, dof_map, constraint_matrix, soln)
    scale_columns(unknowns, input.bies, input.params)
    return unknowns

def dense_matrix(tbem, input, dof_map, constraint_matrix, systems):
    matrix = uncondensed_dense_matrix(tbem, input, dof_map, constraint_matrix, systems)
    condensed_op = tbem.condense_matrix(constraint_matrix, constraint_matrix, matrix)
    op_data = condensed_op.data()
    n_condensed = np.sqrt(op_data.size)
    np_op = op_data.reshape((n_condensed, n_condensed))
    return np_op

def uncondensed_dense_matrix(tbem, input, dof_map, constraint_matrix, systems):
    n = dof_map['n_total_dofs']
    matrix = np.empty((n, n))

    scalings = dict()
    for bie in input.bies:
        scalings[bie['obs_mesh']] = bie['scaling'](input.params)

    for s, bie in zip(systems, input.bies):
        for op in s['lhs']:
            row_components = dof_map[(op['spec']['obs_mesh'], bie['unknown_field'])]
            col_components = dof_map[(op['spec']['src_mesh'], op['spec']['function'])]
            start_row = row_components[0]
            end_row = row_components[-1]
            start_col = col_components[0]
            end_col = col_components[-1]
            rows = end_row - start_row
            cols = end_col - start_col
            reshaped_op = op['op'].data().reshape((rows, cols))
            row_scaling = scalings[bie['obs_mesh']]
            col_scaling = scalings[op['spec']['src_mesh']]
            reshaped_op *= row_scaling * col_scaling
            reshaped_op *= op['spec']['multiplier']
            matrix[start_row:end_row, start_col:end_col] = reshaped_op
    matrix = tbem.DenseOperator(n, n, matrix.reshape(n * n))
    return matrix

def dense_rhs(tbem, input, dof_map, constraint_matrix, systems):
    right_hand_sides = uncondensed_dense_rhs(
        tbem, input, dof_map, constraint_matrix, systems
    )
    rhs = concatenate_condense(tbem, dof_map, constraint_matrix, right_hand_sides)
    return rhs

def uncondensed_dense_rhs(tbem, input, dof_map, constraint_matrix, systems):
    right_hand_sides = [s['rhs'] for s in systems]
    scale_rows(right_hand_sides, input.bies, input.params)
    return right_hand_sides

