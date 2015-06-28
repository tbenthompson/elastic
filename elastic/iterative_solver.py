from dof_handling import distribute_expand, concatenate_condense, \
    scale_columns, scale_rows, operate_on_solution_fields
import numpy as np
import scipy.sparse.linalg

def build_matrix_vector_product(tbem, input, dof_map, constraint_matrix, systems):
    def mat_vec(v):
        unknowns = distribute_expand(tbem, dof_map, constraint_matrix, v)
        scale_columns(unknowns, input.bies, input.params)
        eval = operate_on_solution_fields(tbem, systems, unknowns)
        scale_rows(eval, input.bies, input.params)
        out = concatenate_condense(tbem, dof_map, constraint_matrix, eval)

        if mat_vec.inhomogeneous_constraint_component is not None:
            out -= mat_vec.inhomogeneous_constraint_component

        print("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out
    mat_vec.n_its = 0
    mat_vec.inhomogeneous_constraint_component = None
    return mat_vec


def iterative_solver(tbem, input, dof_map, constraint_matrix, systems):
    #TODO: Do an ILU preconditioning
    mat_vec = build_matrix_vector_product(
        tbem, input, dof_map, constraint_matrix, systems
    )

    def residual_callback(resid):
        print("residual: " + str(resid))
        pass

    right_hand_sides = [s['rhs'] for s in systems]
    scale_rows(right_hand_sides, input.bies, input.params)
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
        tol = input.params['solver_tol'],
        callback = residual_callback
    )
    assert(res[1] == 0)
    unknowns = distribute_expand(tbem, dof_map, constraint_matrix, res[0])
    scale_columns(unknowns, input.bies, input.params)
    return unknowns
