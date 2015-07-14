import numpy as np
from constraints import distribute, condense
from system import scale
import bie_spec
import iterative_solver
import logging
logger = logging.getLogger(__name__)

def dense_solver(tbem, params, dof_map, constraint_matrix, systems):
    log_start_solve(dof_map.n_total_dofs)

    rhs = iterative_solver.calculate_rhs(
        tbem, params, dof_map, constraint_matrix, systems
    )

    homogenized_cm = tbem.homogenize_constraints(constraint_matrix)
    np_op = dense_matrix(tbem, params, dof_map, homogenized_cm, systems)
    log_condition_number(np_op)
    # plot_matrix(np_op)
    soln = np.linalg.solve(np_op, rhs)

    unknowns = iterative_solver.handle_solution(
        tbem, homogenized_cm, dof_map, params, soln
    )
    iterative_solver.add_bcs(tbem, constraint_matrix, dof_map, unknowns)

    log_finish_solve()
    return unknowns

def log_start_solve(n_dofs):
    logger.info('Dense linear solve for system with ' + str(n_dofs) + ' rows.')

def log_condition_number(matrix):
    cond = np.linalg.cond(matrix)
    logger.info('Linear system has condition number: ' + str(cond))
    if cond > 10e12:
        logger.warning('Linear system has very large condition number: ' + str(cond))

def log_finish_solve():
    logger.info('Finished dense linear solve')

def dense_matrix(tbem, params, dof_map, constraint_matrix, systems):
    matrix = uncondensed_dense_matrix(tbem, params, dof_map, systems)
    condensed_op = tbem.condense_matrix(constraint_matrix, constraint_matrix, matrix)
    op_data = condensed_op.data()
    n_condensed = np.sqrt(op_data.size)
    np_op = op_data.reshape((n_condensed, n_condensed))
    return np_op

def uncondensed_dense_matrix(tbem, params, dof_map, systems):
    n = dof_map.n_total_dofs
    matrix = np.zeros((n, n))

    for s in systems:
        for op in s.terms:
            if op.input_type()[1] == 'ones':
                continue
            start_row, end_row, start_col, end_col = dof_map.get_matrix_block(
                s.output_type(), op.input_type()
            )
            row_scaling = bie_spec.integral_scaling[s.output_type()](params)
            col_scaling = bie_spec.integral_scaling[op.input_type()](params)
            chunk = op.to_numpy_matrix() * row_scaling * col_scaling
            log_dense_system_insertion(
                start_row, end_row, start_col, end_col, s.output_type(),
                op.spec, row_scaling * col_scaling
            )
            matrix[start_row:end_row, start_col:end_col] += chunk
    matrix = tbem.DenseOperator(n, n, matrix.reshape(n * n))
    return matrix

def log_skipped_mass(output_type, mass_spec):
    logger.debug('Skipping mass operator: ' + str((output_type, mass_spec)))

def log_dense_system_insertion(*data):
    logger.debug('Insert operator into dense system: ' + str(data))

def plot_matrix(np_op):
    import matplotlib.pyplot as plt
    plt.figure(facecolor = (0, 0, 0))
    plt.imshow(np.log10(np.abs(np_op)), vmin = 3.0, vmax = 8.5)
    plt.axis('off')
    plt.gca().set_axis_bgcolor((0, 0, 0))
    plt.savefig('matrix.pdf', facecolor = (0, 0, 0))

