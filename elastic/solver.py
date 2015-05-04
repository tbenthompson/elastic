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
    for bie in bies:
        obs_mesh = meshes[bie['obs_mesh']]
        for d in range(tbem.dim):
            components.append(obs_mesh.n_dofs())
    return tbempy.build_block_dof_map(components)

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
    return tbempy.from_constraints(out)

def concatenate_condense(dof_map, constraint_matrix, fncs):
    flattened = []
    for f in fncs:
        for d in range(f.size()):
            flattened.append(f.storage[d])
    concat_result = tbempy.concatenate(dof_map, tbempy.BlockVectorX(flattened))
    return tbempy.condense_vector(constraint_matrix, concat_result)

def solve(tbem, input, dof_map, constraint_matrix, systems):
    #TODO: Do a jacobi preconditioning?

    def mat_vec(v):
        vec_v = tbempy.VectorX(v)
        distributed = tbempy.distribute_vector(constraint_matrix, vec_v, dof_map.n_dofs)
        unstacked_soln = tbempy.expand(dof_map, distributed)
        unknowns = extract_solution_fields(tbem, unstacked_soln, input.bies)
        eval = operate_on_solution_fields(tbem, systems, unknowns)
        out = concatenate_condense(dof_map, constraint_matrix, eval)

        print("iteration: " + str(mat_vec.n_its))
        mat_vec.n_its += 1
        return out.storage
    mat_vec.n_its = 0

    def residual_callback(resid):
        print("residual: " + str(resid))

    right_hand_sides = [s['rhs'] for s in systems]
    rhs = concatenate_condense(dof_map, constraint_matrix, right_hand_sides)
    np_rhs = rhs.storage

    A = scipy.sparse.linalg.LinearOperator(
        (np_rhs.shape[0], np_rhs.shape[0]),
        matvec = mat_vec,
        dtype = np.float64
    )

    res = scipy.sparse.linalg.gmres(
        A,
        np_rhs,
        tol = input.params['solver_tol'],
        callback = residual_callback
    )
    assert(res[1] == 0)
    soln = tbempy.VectorX(res[0])
    distributed_soln = tbempy.distribute_vector(constraint_matrix, soln, dof_map.n_dofs)
    unstacked_soln = tbempy.expand(dof_map, distributed_soln)
    unknowns = extract_solution_fields(tbem, unstacked_soln, input.bies)
    return unknowns

def extract_solution_fields(tbem, concatenated_fields, bies):
    fields = dict()
    for i in range(len(bies)):
        vec = []
        for d in range(tbem.dim):
            vec.append(concatenated_fields.storage[i * tbem.dim + d])
        # The field is defined over the obs_mesh of the boundary integral equation
        mesh = bies[i]['obs_mesh']
        name = bies[i]['unknown_field']
        fields[(mesh, name)] = tbempy.BlockVectorX(vec)
    return fields

def operate_on_solution_fields(tbem, systems, unknowns):
    eval = []
    for s in systems:
        evaluated = compute.evaluate_computable_terms(tbem, s['lhs'], unknowns)
        assert(len(evaluated['lhs']) == 0)
        eval.append(-evaluated['rhs'])
    return eval

if __name__ == "__main__":
    Element = input_builder.Element
    soln = controller(2, [
        Element(
            pts = [[0, 0], [1, 0]],
            bc = [[0, 0], [0, 0]],
            bc_type = 'displacement',
            n_refines = 0
        ),
        Element(
            pts = [[1, 0], [2, 0]],
            bc = [[0.005, 0], [0.005, 0]],
            bc_type = 'displacement',
            n_refines = 0
        )
    ], dict(
    ))
    print(soln[('displacement', 'traction')].storage[0].storage)


