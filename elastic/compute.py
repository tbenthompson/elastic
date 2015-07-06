import numpy as np

def form_linear_systems(tbem, bies, meshes, bcs, kernels, params):
    systems = []
    bc_fields = fields_from_bcs(bcs)
    for b in bies:
        setup_bie = setup_integral_equation(tbem, meshes, kernels, params, b)
        s = evaluate_computable_terms(tbem, setup_bie, bc_fields)
        systems.append(s)
    return systems

def setup_integral_equation(tbem, meshes, kernels, params, bie):
    integrals = []
    for term in bie['terms']:
        integrals.append(compute_boundary_integral(tbem, meshes, kernels, params, term))
    integrals.append(compute_mass(tbem, meshes, params, bie['mass_term']))
    return integrals

def make_integrator(tbem, params, kernel):
    return tbem.make_sinh_integrator(
        params['sinh_order'], params['obs_order'],
        params['singular_steps'], params['far_threshold'],
        kernel
    )

def compute_boundary_integral(tbem, meshes, kernels, params, op_spec):
    obs_mesh = meshes[op_spec['obs_mesh']]
    src_mesh = meshes[op_spec['src_mesh']]
    kernel = kernels[op_spec['kernel']]

    mthd = make_integrator(tbem, params, kernel)

    def fmm_boundary_operator(obs_mesh, src_mesh, mthd, all_mesh):
        fmm_config = tbem.FMMConfig(0.3, params['fmm_order'], 250, 0.1, True)
        return tbem.boundary_operator(obs_mesh, src_mesh, mthd, fmm_config, all_mesh)
    op_builder = fmm_boundary_operator
    if params['dense']:
        op_builder = tbem.dense_boundary_operator
    op = op_builder(obs_mesh, src_mesh, mthd, meshes['all_mesh'])
    return dict(op = op, spec = op_spec)

def compute_interior_integral(tbem, meshes, kernels, params,
        op_spec, pts, normals, fields):
    src_mesh = meshes[op_spec['src_mesh']]
    kernel = kernels[op_spec['kernel']]

    mthd = make_integrator(tbem, params, kernel)
    op = tbem.dense_interior_operator(
        pts, normals, src_mesh, mthd, meshes['all_mesh']
    )

    return dict(op = op, spec = op_spec)

def apply_op(op, fields):
    f = fields.get((op['spec']['src_mesh'], op['spec']['function']), None)
    if f is None:
        return None
    computed = op['op'].apply(np.concatenate(f))
    # Negate to switch the term from the LHS to the RHS
    return -computed * op['spec']['multiplier']

def compute_mass(tbem, meshes, params, mass_spec):
    obs_mesh = meshes[mass_spec['obs_mesh']]
    op = tbem.mass_operator_tensor(obs_mesh, params['obs_order'])
    return dict(op = op, spec = dict(
        src_mesh = mass_spec['obs_mesh'],
        function = mass_spec['function'],
        multiplier = mass_spec['multiplier']
    ))

def fields_from_bcs(bcs):
    fields = dict()
    for k, v in bcs.iteritems():
        n_dofs = v.shape[0] * v.shape[1]
        n_components = v.shape[2]
        reshaped = v.reshape((n_dofs, n_components))
        f = [0] * n_components
        for d in range(n_components):
            f[d] = reshaped[:, d]
        fields[(k, k)] = f
    return fields

def evaluate_computable_terms(tbem, bie, fields):
    n_out_dofs = bie[0]['op'].n_rows()
    n_components = tbem.dim
    n_dofs_per_component = n_out_dofs / n_components;
    evaluated_terms = np.zeros(tbem.dim * n_dofs_per_component)
    uncomputed_terms = []
    for t in bie:
        computed = apply_op(t, fields)
        if computed is not None:
            evaluated_terms += computed
        else:
            uncomputed_terms.append(t)
    return dict(lhs = uncomputed_terms, rhs = evaluated_terms)

def interior_eval(tbem, bcs, meshes, kernels, params, soln, pts, normals, terms):
    for k, v in soln.iteritems():
        fields[k] = v
    result = np.zeros(pts.shape[0] * tbem.dim)
    for t in terms:
        op = compute_interior_integral(
            tbem, meshes, kernels, params, t, pts, normals, fields
        )
        result += apply_op(op, fields)
    out = []
    for d in range(tbem.dim):
        start_idx = d * pts.shape[0]
        end_idx = (d + 1) * pts.shape[0]
        out.append(result[start_idx:end_idx])
    return out
