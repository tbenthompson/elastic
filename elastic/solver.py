import scipy.sparse.linalg
import numpy as np

import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import input_builder
import compute


class Result(object):
    def __init__(self, tbem, soln, input):
        self.tbem = tbem
        self.soln = soln
        self.input = input

    def interior_displacement(self, pts):
        bie = bie_spec.get_displacement_BIE("displacement", self.input.params)
        normals = np.zeros_like(pts)
        return self.interior_eval(pts, normals, bie)

    def interior_traction(self, pts, normals):
        bie = bie_spec.get_traction_BIE("traction", self.input.params)
        return self.interior_eval(pts, normals, bie)

    def interior_eval(self, pts, normals, bie):
        fields = compute.fields_from_bcs(self.input.bcs)
        for k, v in self.soln.iteritems():
            fields[k] = v
        result = np.zeros(pts.shape[0] * self.tbem.dim)
        for term in bie['terms']:
            src_mesh = self.input.meshes[term['src_mesh']]
            kernel = self.input.kernels[term['kernel']]
            f = fields[(term['src_mesh'], term['function'])]

            mthd = self.tbem.make_sinh_integrator(
                self.input.params['sinh_order'], self.input.params['obs_order'],
                self.input.params['singular_steps'], self.input.params['far_threshold'],
                kernel
            )
            op = self.tbem.dense_interior_operator(
                pts, normals, src_mesh, mthd, self.input.all_mesh
            )

            # The result is negated to move it to the other side of the equation
            result -= op.apply(np.concatenate(f)) * term['multiplier']
        out = []
        for d in range(self.tbem.dim):
            start_idx = d * pts.shape[0]
            end_idx = (d + 1) * pts.shape[0]
            out.append(result[start_idx:end_idx])
        return out

    def save(self, filename):
        with open(filename, 'w') as f:
            np.savez(
                f,
                dim = self.tbem.dim,
                soln = self.soln,
                elements = self.input.elements,
                params = self.input.params
            )

    @staticmethod
    def load(filename):
        with open(filename, 'r') as f:
            npzfile = np.load(f)
            soln = npzfile['soln'].tolist()
            raw_elements = npzfile['elements'].tolist()
            elements = [
                input_builder.Element(*r) for r in raw_elements
            ]
            params = npzfile['params'].tolist()
            dim = npzfile['dim']
        tbem = get_tbem(dim)
        input = input_builder.build_input(tbem, elements, params)
        return Result(tbem, soln, input)

def execute(dim, elements, input_params):
    tbem, input, dof_map, constraint_matrix, systems = \
        form_system(dim, elements, input_params)
    solve_fnc = iterative_solver
    if input.params['dense']:
        solve_fnc = dense_solver
    soln = solve_fnc(
        tbem, input, dof_map, constraint_matrix, systems
    )
    return Result(tbem, soln, input)

def form_system(dim, elements, input_params):
    dim = dim
    elements = elements
    input_params = input_params
    tbem = get_tbem(dim)
    input = input_builder.build_input(tbem, elements, input_params)
    dof_map = build_dof_map(tbem, input.bies, input.meshes)
    constraint_matrix = build_constraint_matrix(
        tbem, dof_map, input
    )
    systems = compute.form_linear_systems(tbem, input)
    return tbem, input, dof_map, constraint_matrix, systems

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

def build_constraint_matrix(tbem, dof_map, input):
    out = []
    for bie in input.bies:
        constraints = bie['constraint_builder'](tbem, dof_map, input.meshes, input)
        out.extend(constraints)
    return tbem.from_constraints(out)

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

def dense_matrix(tbem, input, dof_map, constraint_matrix, systems):
    matrix = uncondensed_dense_matrix(tbem, input, dof_map, constraint_matrix, systems)
    condensed_op = tbem.condense_matrix(constraint_matrix, constraint_matrix, matrix)
    op_data = condensed_op.data()
    n_condensed = np.sqrt(op_data.size)
    np_op = op_data.reshape((n_condensed, n_condensed))
    return np_op

def uncondensed_dense_rhs(tbem, input, dof_map, constraint_matrix, systems):
    right_hand_sides = [s['rhs'] for s in systems]
    scale_rows(right_hand_sides, input.bies, input.params)
    return right_hand_sides

def dense_rhs(tbem, input, dof_map, constraint_matrix, systems):
    right_hand_sides = uncondensed_dense_rhs(
        tbem, input, dof_map, constraint_matrix, systems
    )
    rhs = concatenate_condense(tbem, dof_map, constraint_matrix, right_hand_sides)
    return rhs

def dense_solver(tbem, input, dof_map, constraint_matrix, systems):
    np_op = dense_matrix(tbem, input, dof_map, constraint_matrix, systems)
    rhs = dense_rhs(tbem, input, dof_map, constraint_matrix, systems)
    input.logger.linear_solve_start(rhs.shape[0])
    soln = np.linalg.solve(np_op, rhs)
    input.logger.linear_solve_end()
    unknowns = distribute_expand(tbem, dof_map, constraint_matrix, soln)
    scale_columns(unknowns, input.bies, input.params)
    return unknowns

def build_matrix_vector_product(tbem, input, dof_map, constraint_matrix, systems):
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

    print(mat_vec(np.zeros_like(rhs)))

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
