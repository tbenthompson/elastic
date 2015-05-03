from tbempy import BlockVectorX, VectorX

def form_linear_systems(tbem, input):
    systems = []
    for bie in input.bies:
        setup_bie = setup_integral_equation(tbem, input, bie)
        bc_fields = fields_from_bcs(input.bcs)
        s = evaluate_computable_terms(tbem, setup_bie, bc_fields)
        systems.append(s)
    return systems

def setup_integral_equation(tbem, input, bie):
    integrals = []
    for term in bie['terms']:
        integrals.append(compute_integral(tbem, input, term))
    integrals.append(compute_mass(tbem, input, bie['mass_term']))
    return integrals

def compute_integral(tbem, input, op_spec):
    obs_mesh = input.meshes[op_spec['obs_mesh']]
    src_mesh = input.meshes[op_spec['src_mesh']]
    kernel = input.kernels[op_spec['kernel']]

    mthd = tbem.make_adaptive_integration_mthd(input.quad_strategy, kernel)
    op = tbem.integral_operator(obs_mesh, src_mesh, mthd)
    return dict(op = op, spec = op_spec)

def compute_mass(tbem, input, mass_spec):
    obs_mesh = input.meshes[mass_spec['obs_mesh']]
    op = tbem.mass_operator_tensor(obs_mesh, input.params['obs_order'])
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
        x = BlockVectorX([
            VectorX(reshaped[:, d]) for d in range(n_components)
        ])
        fields[(k, k)] = x
    return fields

def evaluate_computable_terms(tbem, bie, fields):
    n_total_out_dofs = bie[0]['op'].n_total_rows()
    n_components = tbem.dim
    n_dofs_per_component = n_total_out_dofs / n_components;
    evaluated_terms = BlockVectorX(tbem.dim, VectorX(n_dofs_per_component, 0.0))
    uncomputed_terms = []
    for t in bie:
        f = fields.get((t['spec']['src_mesh'], t['spec']['function']), None)
        if f is None:
            uncomputed_terms.append(t)
        else:
            #TODO: better numpy to vectorx translation would make this unnecessary
            computed = t['op'].apply(f)
            # Negate the computed value because it is being moved from the LHS
            # to the RHS of the system
            evaluated_terms -= computed * t['spec']['multiplier']
    return dict(lhs = uncomputed_terms, rhs = evaluated_terms)
