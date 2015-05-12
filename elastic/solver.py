import scipy.sparse.linalg
import numpy as np

import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import input_builder
import compute

class execute(object):
    def __init__(self, dim, elements, input_params):
        self.dim = dim
        self.elements = elements
        self.input_params = input_params
        tbem = get_tbem(dim)
        self.input = input_builder.build_input(tbem, self.elements, self.input_params)
        self.dof_map = build_dof_map(tbem, self.input.bies, self.input.meshes)
        self.constraint_matrix = build_constraint_matrix(
            tbem, self.dof_map, self.input.bies, self.input.meshes
        )
        self.systems = compute.form_linear_systems(tbem, self.input)
        solve_fnc = iterative_solver
        if self.input.params['dense']:
            solve_fnc = dense_solver
        self.soln = solve_fnc(
            tbem, self.input, self.dof_map, self.constraint_matrix, self.systems
        )

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def build_dof_map(tbem, bies, meshes):
    dof_map = dict()
    start = 0
    for bie in bies:
        components = []
        n_dofs = meshes[bie['obs_mesh']].n_dofs()
        for d in range(tbem.dim):
            components.append(start)
            start += n_dofs
        components.append(start)
        dof_map[(bie['obs_mesh'], bie['unknown_field'])] = components
    dof_map['n_total_dofs'] = start
    return dof_map

def build_constraint_matrix(tbem, dof_map, bies, meshes):
    out = []
    for bie in bies:
        constraints = bie['constraint_builder'](tbem, dof_map, meshes)
        out.extend(constraints)
    return tbem.from_constraints(out)

def dense_solver(tbem, input, dof_map, constraint_matrix, systems):
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
            matrix[start_row:end_row, start_col:end_col] = reshaped_op
    matrix = tbem.DenseOperator(n, n, matrix.reshape(n * n))
    condensed_op = tbem.condense_matrix(constraint_matrix, constraint_matrix, matrix)
    right_hand_sides = [s['rhs'] for s in systems]
    # TODO: Figure out why this negation is necessary for the dense solver but
    # not the iterative solver!
    right_hand_sides[0] = -right_hand_sides[0]
    scale_rows(right_hand_sides, input.bies, input.params)
    rhs = concatenate_condense(tbem, dof_map, constraint_matrix, right_hand_sides)
    np_op = condensed_op.data().reshape((rhs.shape[0], rhs.shape[0]))
    soln = np.linalg.solve(np_op, rhs)
    unknowns = distribute_expand(tbem, dof_map, constraint_matrix, soln)
    scale_columns(unknowns, input.bies, input.params)
    return unknowns

def iterative_solver(tbem, input, dof_map, constraint_matrix, systems):
    #TODO: Do an ILU preconditioning

    def mat_vec(v):
        unknowns = distribute_expand(tbem, dof_map, constraint_matrix, v)
        scale_columns(unknowns, input.bies, input.params)
        eval = operate_on_solution_fields(tbem, systems, unknowns)
        scale_rows(eval, input.bies, input.params)
        out = concatenate_condense(tbem, dof_map, constraint_matrix, eval)

        print("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out
    mat_vec.n_its = 0

    def residual_callback(resid):
        print("residual: " + str(resid))

    right_hand_sides = [s['rhs'] for s in systems]
    scale_rows(right_hand_sides, input.bies, input.params)
    rhs = concatenate_condense(tbem, dof_map, constraint_matrix, right_hand_sides)

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

def distribute_expand(tbem, dof_map, constraint_matrix, vec):
    distributed = tbem.distribute_vector(constraint_matrix, vec, dof_map['n_total_dofs'])
    assert(distributed.shape[0] == dof_map['n_total_dofs'])
    fields = dict()
    for mesh_and_field,dofs in dof_map.iteritems():
        if mesh_and_field == 'n_total_dofs':
            continue
        field_components = []
        for d in range(tbem.dim):
            start_idx = dofs[d]
            end_idx = dofs[d + 1]
            field_components.append(distributed[start_idx:end_idx])
        fields[mesh_and_field] = field_components
    return fields

def concatenate_condense(tbem, dof_map, constraint_matrix, fncs):
    concat = np.concatenate(fncs)
    return tbem.condense_vector(constraint_matrix, concat)

def operate_on_solution_fields(tbem, systems, unknowns):
    eval = []
    for s in systems:
        evaluated = compute.evaluate_computable_terms(tbem, s['lhs'], unknowns)
        assert(len(evaluated['lhs']) == 0)
        eval.append(-evaluated['rhs'])
    return eval

def scale_columns(unknowns, bies, params):
    for bie in bies:
        o = bie['obs_mesh']
        f = bie['unknown_field']
        values = unknowns[(o, f)]
        for d in range(len(values)):
            values[d] *= bie['scaling'](params)

def scale_rows(eval, bies, params):
    for i, bie in enumerate(bies):
        eval[i] *= bie['scaling'](params)
