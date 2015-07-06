from dof_handling import distribute_expand, concatenate_condense, \
    scale_columns, scale_rows, operate_on_solution_fields

import numpy as np
import scipy.sparse.linalg
import logging
logger = logging.getLogger(__name__)

def build_matrix_vector_product(tbem, params, bies, dof_map, constraint_matrix, systems):
    def mat_vec(v):
        unknowns = distribute_expand(tbem, dof_map, constraint_matrix, v)
        scale_columns(unknowns, bies, params)
        eval = operate_on_solution_fields(tbem, systems, unknowns)
        scale_rows(eval, bies, params)
        out = concatenate_condense(tbem, dof_map, constraint_matrix, eval)

        if mat_vec.inhomogeneous_constraint_component is not None:
            out -= mat_vec.inhomogeneous_constraint_component

        logger.info("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out
    mat_vec.n_its = 0
    mat_vec.inhomogeneous_constraint_component = None
    return mat_vec


def iterative_solver(tbem, params, bies, dof_map, constraint_matrix, systems):
    logger.debug('Iteratively solving linear system')
    #TODO: Do an ILU preconditioning
    mat_vec = build_matrix_vector_product(
        tbem, params, bies, dof_map, constraint_matrix, systems
    )

    def residual_callback(resid):
        logger.info("residual: " + str(resid))
        pass

    right_hand_sides = [s['rhs'] for s in systems]
    scale_rows(right_hand_sides, bies, params)
    rhs = concatenate_condense(tbem, dof_map, constraint_matrix, right_hand_sides)

    inhomogeneous = mat_vec(np.zeros_like(rhs))
    rhs -= inhomogeneous
    mat_vec.inhomogeneous_constraint_component = inhomogeneous

    A = scipy.sparse.linalg.LinearOperator(
        (rhs.shape[0], rhs.shape[0]),
        matvec = mat_vec,
        dtype = np.float64
    )

    res = scipy.sparse.linalg.gmres(
        A,
        rhs,
        tol = params['solver_tol'],
        callback = residual_callback
    )
    assert(res[1] == 0)
    unknowns = distribute_expand(tbem, dof_map, constraint_matrix, res[0])
    scale_columns(unknowns, bies, params)
    logger.debug(
        'Iterative solution complete, requiring %d iterations', mat_vec.n_its
    )
    return unknowns
