from collections import namedtuple
import numpy as np

def form_traction_constraints(tbem, dof_map, meshes):
    return []

def form_slip_constraints(tbem, dof_map, meshes):
    return []
    # continuity = tbem.mesh_continuity(meshes['discontinuous'].begin())
    # one_component = tbem.convert_to_constraints(continuity)
    # all_components = []
    # for d in range(tbem.dim):
    #     start_dof = dof_map.get('discontinuous', 'slip', d)
    #     all_components.extend(tbem.shift_constraints(one_component, start_dof))
    # return all_components

def form_displacement_constraints(tbem, dof_map, meshes):
    continuity = tbem.mesh_continuity(meshes['continuous'].begin())
    cut_continuity = tbem.cut_at_intersection(
        continuity,
        meshes['continuous'].begin(),
        meshes['discontinuous'].begin()
    )
    one_component = tbem.convert_to_constraints(cut_continuity)
    all_components = []
    for d in range(tbem.dim):
        start_dof = dof_map.get('continuous', 'displacement', d)
        all_components.extend(tbem.shift_constraints(one_component, start_dof))
    return all_components

def gather_continuity_constraints(tbem, dof_map, meshes):
    return form_traction_constraints(tbem, dof_map, meshes) +\
        form_slip_constraints(tbem, dof_map, meshes) +\
        form_displacement_constraints(tbem, dof_map, meshes)

Term = namedtuple('Term', ['field', 'component', 'weight'])
Constraint = namedtuple('Constraint', ['terms', 'rhs'])
BCConstraint = namedtuple('BCConstraint', ['field', 'component', 'rhs'])

#TODO: Much of this BC constraint handling code would be better in C++...
def gather_bc_constraints(tbem, dof_map, es):
    tbempy_constraints = []
    counts = dict(discontinuous = 0, continuous = 0)
    for e in es:
        #TODO: Handle the bad input instead of failing
        assert(len(e['constraints']) == 2)
        for basis_idx in range(tbem.dim):
            for c in e['constraints']:
                tbempy_constraints.append(convert_helper(
                    tbem, dof_map, counts[e['type']], basis_idx, c, e)
                )
        counts[e['type']] += 1
    return tbempy_constraints

def convert_helper(tbem, dof_map, element_idx, basis_idx, c, element):
    if type(c) is BCConstraint:
        c = convert_bcconstraint_to_constraint(c)

    new_terms = add_element_local_terms(tbem.dim, c.terms, element['pts'])

    tbempy_terms = []
    for t in new_terms:
        component_start = dof_map.get(element['type'], t.field, t.component)
        dof = component_start + tbem.dim * element_idx + basis_idx
        tbempy_terms.append(tbem.LinearTerm(dof, t.weight))

    return tbem.ConstraintEQ(tbempy_terms, c.rhs[basis_idx])

def add_element_local_terms(dim, terms, pts):
    new_terms = []
    for t in terms:
        if type(t.component) is not str:
            new_terms.append(t)
            continue
        new_terms.extend(transform_element_local_term(dim, t, pts))
    return new_terms

def transform_element_local_term(dim, term, pts):
    if term.component == 'tangential0':
        direction = tangential_vector0(pts)
    elif term.component == 'tangential1' and tbem.dim == 3:
        return 'unimplemented!'
    elif term.component == 'normal':
        direction = normal_vector(pts)
    else:
        return 'Not a valid element local direction.'
    return [Term(term.field, d, direction[d] * term.weight) for d in range(dim)]

def tangential_vector0(pts):
    vector = np.array([
        pts[1][0] - pts[0][0],
        pts[1][1] - pts[0][1]
    ])
    return vector / np.linalg.norm(vector)

def normal_vector(pts):
    vector = np.array([
        -(pts[1][1] - pts[0][1]),
        pts[1][0] - pts[0][0]
    ])
    return vector / np.linalg.norm(vector)

def build_constraint_matrix(tbem, dof_map, elements, meshes):
    bc_constraints = gather_bc_constraints(tbem, dof_map, elements)
    continuity_constraints = gather_continuity_constraints(tbem, dof_map, meshes)
    all = bc_constraints + continuity_constraints
    constraint_matrix = tbem.from_constraints(bc_constraints + continuity_constraints)
    return constraint_matrix

def distribute(tbem, constraint_matrix, n_total_dofs, vec):
    distributed = tbem.distribute_vector(constraint_matrix, vec, n_total_dofs)
    assert(distributed.shape[0] == n_total_dofs)
    return distributed

def condense(tbem, constraint_matrix, concatenated):
    return tbem.condense_vector(constraint_matrix, concatenated)

def convert_bcconstraint_to_constraint(c):
    return Constraint(terms = [Term(c.field, c.component, 1.0)], rhs = c.rhs)

def refine_constraint(c, end_idx):
    if type(c) is BCConstraint:
        c = convert_bcconstraint_to_constraint(c)
    midpt_rhs = (c.rhs[0] + c.rhs[1]) / 2.0
    if end_idx == 0:
        rhs = [c.rhs[0], midpt_rhs]
    else:
        rhs = [midpt_rhs, c.rhs[0]]
    return Constraint(terms = c.terms, rhs = rhs)

