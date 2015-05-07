import scipy.sparse.linalg
import numpy as np

import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import input_builder
import compute

class Controller(object):
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
        self.soln = solve(
            tbem, self.input, self.dof_map, self.constraint_matrix, self.systems
        )

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def build_dof_map(tbem, bies, meshes):
    components = []
    start = 0
    for bie in bies:
        obs_mesh = meshes[bie['obs_mesh']]
        for d in range(tbem.dim):
            components.append(start)
            n_dofs = obs_mesh.n_dofs()
            start += n_dofs
    components.append(start)
    return components

def build_constraint_matrix(tbem, dof_map, bies, meshes):
    eqtn_component_map = dict()
    for i in range(len(bies)):
        eqtn_component_map[bies[i]['obs_mesh']] = i * tbem.dim

    out = []
    for bie in bies:
        for d in range(tbem.dim):
            constraints = bie['constraint_builder'](
                tbem, eqtn_component_map, dof_map, meshes, d
            )
            out.extend(constraints)
    return tbem.from_constraints(out)

def distribute_expand(tbem, dof_map, constraint_matrix, vec):
    distributed = tbem.distribute_vector(constraint_matrix, vec, dof_map[-1])
    assert(distributed.shape[0] == dof_map[-1])
    unstacked_soln = []
    for i in range(len(dof_map) - 1):
        unstacked_soln.append(distributed[dof_map[i]:dof_map[i+1]])
    return unstacked_soln

def concatenate_condense(tbem, dof_map, constraint_matrix, fncs):
    concat = np.concatenate(fncs)
    return tbem.condense_vector(constraint_matrix, concat)

def solve(tbem, input, dof_map, constraint_matrix, systems):
    #TODO: Do a jacobi preconditioning?

    def mat_vec(v):
        unstacked_soln = distribute_expand(tbem, dof_map, constraint_matrix, v)
        unknowns = extract_solution_fields(tbem, unstacked_soln, input.bies)
        eval = operate_on_solution_fields(tbem, systems, unknowns)
        out = concatenate_condense(tbem, dof_map, constraint_matrix, eval)

        print("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out
    mat_vec.n_its = 0

    def residual_callback(resid):
        print("residual: " + str(resid))

    right_hand_sides = [s['rhs'] for s in systems]
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
    unstacked_soln = distribute_expand(tbem, dof_map, constraint_matrix, res[0])
    unknowns = extract_solution_fields(tbem, unstacked_soln, input.bies)
    return unknowns

def extract_solution_fields(tbem, concatenated_fields, bies):
    fields = dict()
    for i in range(len(bies)):
        vec = []
        for d in range(tbem.dim):
            vec.append(concatenated_fields[i * tbem.dim + d])
        # The field is defined over the obs_mesh of the boundary integral equation
        mesh = bies[i]['obs_mesh']
        name = bies[i]['unknown_field']
        fields[(mesh, name)] = vec
    return fields

def operate_on_solution_fields(tbem, systems, unknowns):
    eval = []
    for s in systems:
        evaluated = compute.evaluate_computable_terms(tbem, s['lhs'], unknowns)
        assert(len(evaluated['lhs']) == 0)
        eval.append(-evaluated['rhs'])
    return eval

