import bie_spec
from constraints import distribute, condense
from compute import evaluate_linear_systems, scale

import numpy as np
import scipy.sparse.linalg
import logging
logger = logging.getLogger(__name__)

def info(abc):
    print(np.max(abc[('continuous', 'traction')]))
    print(np.max(abc[('continuous', 'displacement')]))
scale_factor = 4000.0

def build_matrix_vector_product(tbem, params, bies, dof_map, constraint_matrix, systems):
    def mat_vec(v):
        distributed = distribute(tbem, constraint_matrix, dof_map.n_total_dofs, v)
        unknowns = dof_map.expand(distributed)
        # scale(unknowns, bie_spec.solution_scaling, params, False)
        for d in range(2):
            unknowns[('continuous', 'displacement')][d] /= 1.0
            unknowns[('continuous', 'traction')][d] *= scale_factor

        eval = evaluate_linear_systems(systems, unknowns)
        for d in range(2):
            eval[('continuous', 'displacement')][d] /= 1.0
            eval[('continuous', 'traction')][d] /= scale_factor
        # scale(eval, bie_spec.solution_scaling, params, True)
        concatenated = dof_map.concatenate(eval)
        out = condense(tbem, constraint_matrix, concatenated)
        logger.debug("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out
    mat_vec.n_its = 0
    return mat_vec

def calculate_rhs(tbem, params, bies, dof_map, constraint_matrix, systems):
    rhs_mat_vec = build_matrix_vector_product(
        tbem, params, bies, dof_map, constraint_matrix, systems
    )
    return rhs_mat_vec(np.zeros(dof_map.n_total_dofs))

def iterative_solver(tbem, params, bies, dof_map, constraint_matrix, systems):
    logger.info('Iteratively solving linear system')
    #TODO: Do an ILU preconditioning
    rhs = -calculate_rhs(tbem, params, bies, dof_map, constraint_matrix, systems)
    import matplotlib.pyplot as plt
    plt.plot(rhs)
    plt.show()

    homogenized_cm = tbem.homogenize_constraints(constraint_matrix)
    mat_vec = build_matrix_vector_product(
        tbem, params, bies, dof_map, homogenized_cm, systems
    )
    A = scipy.sparse.linalg.LinearOperator(
        (rhs.shape[0], rhs.shape[0]),
        matvec = mat_vec,
        dtype = np.float64
    )
    res = scipy.sparse.linalg.gmres(
        A,
        rhs,
        tol = params['solver_tol'],
        callback = lambda resid: logger.debug("residual: " + str(resid))
    )
    assert(res[1] == 0)

    distributed = distribute(tbem, homogenized_cm, dof_map.n_total_dofs, res[0])
    unknowns = dof_map.expand(distributed)
    for d in range(2):
        unknowns[('continuous', 'displacement')][d] *= 1.0
        # unknowns[('continuous', 'traction')][d] /= scale_factor
    # scale(unknowns, bie_spec.integral_scaling, params, True)

    inhomogeneous_component = distribute(
        tbem, constraint_matrix, dof_map.n_total_dofs,
        np.zeros(dof_map.n_total_dofs)
    )
    inhomogeneous_unknowns = dof_map.expand(inhomogeneous_component)
    import matplotlib.pyplot as plt
    plt.plot(unknowns[('continuous', 'traction')][0], 'k')
    plt.plot(inhomogeneous_unknowns[('continuous', 'traction')][0], 'r')
    plt.show()
    for k in unknowns:
        for d in range(tbem.dim):
            unknowns[k][d] += inhomogeneous_unknowns[k][d]

    logger.info(
        'Iterative solution complete, requiring %d iterations', mat_vec.n_its
    )
    return unknowns
