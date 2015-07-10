from collections import namedtuple

def form_traction_constraints(tbem, dof_map, meshes):
    return []

def form_slip_constraints(tbem, dof_map, meshes):
    return []

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
        c = Constraint(terms = [Term(c.field, c.component, 1.0)], rhs = c.rhs)

    tbempy_terms = []
    for t in c.terms:
        component_start = dof_map.get(element['type'], t.field, t.component)
        dof = component_start + tbem.dim * element_idx + basis_idx
        tbempy_terms.append(tbem.LinearTerm(dof, t.weight))
    assert(len(tbempy_terms) == len(c.terms))

    return tbem.ConstraintEQ(tbempy_terms, c.rhs[basis_idx])

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
