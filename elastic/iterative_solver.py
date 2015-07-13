import bie_spec
import constraints
import system

import numpy as np
import scipy.sparse.linalg
import logging
logger = logging.getLogger(__name__)

def build_matrix_vector_product(
    tbem, params, dof_map, constraint_matrix, systems, compute_rhs = False):
    def mat_vec(v):
        unknowns = dof_map.expand(constraints.distribute(
            tbem, constraint_matrix, dof_map.n_total_dofs, v
        ))
        if not compute_rhs:
            system.scale(unknowns, bie_spec.integral_scaling, params, False)
        system.add_constant_fields(unknowns, not compute_rhs)
        eval = system.evaluate_linear_systems(systems, unknowns)
        system.scale(eval, bie_spec.integral_scaling, params, False)
        concatenated = dof_map.concatenate(eval)
        out = constraints.condense(tbem, constraint_matrix, concatenated)
        logger.debug("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out
    mat_vec.n_its = 0
    return mat_vec

def calculate_rhs(tbem, params, dof_map, constraint_matrix, systems):
    rhs_mat_vec = build_matrix_vector_product(
        tbem, params, dof_map, constraint_matrix, systems, compute_rhs = True
    )
    return -rhs_mat_vec(np.zeros(dof_map.n_total_dofs))

def handle_solution(tbem, homogenized_cm, dof_map, params, soln):
    unknowns = dof_map.expand(constraints.distribute(
        tbem, homogenized_cm, dof_map.n_total_dofs, soln
    ))
    system.scale(unknowns, bie_spec.integral_scaling, params, False)
    return unknowns

def add_bcs(tbem, constraint_matrix, dof_map, unknowns):
    inhomogeneous_component = constraints.distribute(
        tbem, constraint_matrix, dof_map.n_total_dofs,
        np.zeros(dof_map.n_total_dofs)
    )
    inhomogeneous_unknowns = dof_map.expand(inhomogeneous_component)
    for k in unknowns:
        for d in range(tbem.dim):
            unknowns[k][d] += inhomogeneous_unknowns[k][d]

def gmres_wrapper(mat_vec_fnc, rhs, tol):
    A = scipy.sparse.linalg.LinearOperator(
        (rhs.shape[0], rhs.shape[0]),
        matvec = mat_vec_fnc,
        dtype = np.float64
    )
    res = scipy.sparse.linalg.gmres(
        A,
        rhs,
        tol = tol,
        callback = lambda resid: logger.debug("residual: " + str(resid))
    )
    assert(res[1] == 0)
    return res[0]

def iterative_solver(tbem, params, dof_map, constraint_matrix, systems):
    logger.info('Iteratively solving linear system')
    rhs = calculate_rhs(tbem, params, dof_map, constraint_matrix, systems)

    homogenized_cm = tbem.homogenize_constraints(constraint_matrix)
    mat_vec = build_matrix_vector_product(
        tbem, params, dof_map, homogenized_cm, systems
    )
    soln = gmres_wrapper(mat_vec, rhs, params['solver_tol'])
    unknowns = handle_solution(tbem, homogenized_cm, dof_map, params, soln)
    add_bcs(tbem, constraint_matrix, dof_map, unknowns)

    logger.info(
        'Iterative solution complete, requiring %d iterations', mat_vec.n_its
    )
    return unknowns
