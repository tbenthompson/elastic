import numpy as np

import compute

def build_constraint_matrix(tbem, dof_map, input):
    out = []
    for bie in input.bies:
        constraints = bie['constraint_builder'](tbem, dof_map, input.meshes, input)
        out.extend(constraints)
    return tbem.from_constraints(out)

def distribute_expand(tbem, dof_map, constraint_matrix, vec):
    distributed = tbem.distribute_vector(constraint_matrix, vec, dof_map['n_total_dofs'])
    assert(distributed.shape[0] == dof_map['n_total_dofs'])
    fields = dict()
    for mesh_and_field, dofs in dof_map.iteritems():
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
