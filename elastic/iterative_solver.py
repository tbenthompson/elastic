import bie_spec
from constraints import distribute, condense
from dof_handling import scale_columns, scale_rows
from compute import evaluate_linear_systems

import numpy as np
import scipy.sparse.linalg
import logging
logger = logging.getLogger(__name__)

def build_matrix_vector_product(tbem, params, bies, dof_map, constraint_matrix, systems):
    def mat_vec(v):
        distributed = distribute(tbem, constraint_matrix, dof_map.n_total_dofs, v)
        unknowns = dof_map.expand(distributed)
        scale_columns(unknowns, bie_spec.field_types, params)
        eval = evaluate_linear_systems(systems, unknowns)
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

    inhomogeneous = mat_vec(np.zeros(dof_map.n_total_dofs))
    rhs =- inhomogeneous
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
