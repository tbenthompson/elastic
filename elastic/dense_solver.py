import numpy as np
from constraints import distribute, condense
from system import scale
import bie_spec
import iterative_solver
from log_tools import log_elapsed_time
import logging
logger = logging.getLogger(__name__)

class DenseSolver(object):
    def __init__(self, tbem, params):
        self.tbem = tbem
        self.params = params

    def solve(self, dof_map, constraint_matrix, systems):
        rhs = iterative_solver.calculate_rhs(
            self.tbem, self.params, dof_map, constraint_matrix, systems
        )

        homogenized_cm = self.tbem.homogenize_constraints(constraint_matrix)
        uncondensed_matrix = self.form_dense_matrix(dof_map, systems)
        matrix = self.condense_matrix(uncondensed_matrix, constraint_matrix)
        soln = self.numpy_solve(matrix, rhs)

        unknowns = iterative_solver.handle_solution(
            self.tbem, homogenized_cm, dof_map, self.params, soln
        )
        iterative_solver.add_bcs(self.tbem, constraint_matrix, dof_map, unknowns)
        return unknowns

    @log_elapsed_time(logger, 'insertion of operators into dense system')
    def form_dense_matrix(self, dof_map, systems):
        ops = []
        start_rows = []
        start_cols = []
        multipliers = []
        for s in systems:
            for op in s.terms:
                if op.input_type()[1] == 'ones':
                    continue
                start_row, end_row, start_col, end_col = dof_map.get_matrix_block(
                    s.output_type(), op.input_type()
                )
                row_scaling = bie_spec.integral_scaling[s.output_type()](self.params)
                col_scaling = bie_spec.integral_scaling[op.input_type()](self.params)
                mult = op.get_multiplier() * row_scaling * col_scaling

                ops.append(op.to_dense())
                start_rows.append(start_row)
                start_cols.append(start_col)
                multipliers.append(mult)
                self.log_dense_system_insertion(
                    start_row, end_row, start_col, end_col, s.output_type(),
                    op.spec, row_scaling * col_scaling
                )

        n = dof_map.n_total_dofs
        matrix = self.tbem.compose_dense_ops(
            ops, start_rows, start_cols, multipliers, n, n
        )
        return matrix

    def log_skipped_mass(self, output_type, mass_spec):
        logger.debug('Skipping mass operator: ' + str((output_type, mass_spec)))

    def log_dense_system_insertion(self, *data):
        logger.debug('Insert operator into dense system: ' + str(data))

    @log_elapsed_time(logger, 'condensation of dense system')
    def condense_matrix(self, uncondensed_matrix, constraint_matrix):
        condensed_op = self.tbem.condense_matrix(
            constraint_matrix, constraint_matrix, uncondensed_matrix
        )
        op_data = condensed_op.data()
        n_condensed = np.sqrt(op_data.size)
        matrix = op_data.reshape((n_condensed, n_condensed))
        return matrix

    def log_condition_number_if_desired(self, matrix):
        if not self.params['check_condition_number']:
            return
        cond = np.linalg.cond(matrix)
        logger.info('Linear system has condition number: ' + str(cond))
        if cond > 10e12:
            logger.warning(
                'Linear system has very large condition number: ' + str(cond)
            )

    def log_constrained_system_size(self, matrix):
        logger.info('Constrained linear system has ' + str(matrix.shape[0]) + ' rows.')

    @log_elapsed_time(logger, 'numpy linear solve')
    def numpy_solve(self, matrix, rhs):
        self.log_constrained_system_size(matrix)
        self.log_condition_number_if_desired(matrix)
        # plot_matrix(matrix)
        return np.linalg.solve(matrix, rhs)



def plot_matrix(np_op):
    import matplotlib.pyplot as plt
    plt.figure(facecolor = (0, 0, 0))
    plt.imshow(np.log10(np.abs(np_op)), vmin = 3.0, vmax = 8.5)
    plt.axis('off')
    plt.gca().set_axis_bgcolor((0, 0, 0))
    plt.savefig('matrix.pdf', facecolor = (0, 0, 0))
