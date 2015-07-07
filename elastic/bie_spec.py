
mesh_types = [
    'continuous',
    'discontinuous'
]

unknowns_to_knowns = dict(
    traction = 'displacement',
    displacement = 'traction',
    crack_traction = 'slip',
    slip = 'crack_traction'
)

field_to_mesh = dict(
    displacement = 'continuous',
    traction = 'continuous',
    slip = 'discontinuous',
    crack_traction = 'discontinuous'
)

'''
Returns the elastic kernels (also called Green's functions or
fundamental solutions) for the Somigliana identity and the hypersingular
integral equation
'''
def get_elastic_kernels(tbem, params):
    mu = params['shear_modulus']
    pr = params['poisson_ratio']
    g = params['gravity_vector']
    kernels = dict(
        displacement = tbem.ElasticDisplacement(mu, pr),
        traction = tbem.ElasticTraction(mu, pr),
        adjoint_traction = tbem.ElasticAdjointTraction(mu, pr),
        hypersingular = tbem.ElasticHypersingular(mu, pr),
        gravity_displacement = tbem.GravityDisplacement(mu, pr, g),
        gravity_traction = tbem.GravityTraction(mu, pr, g)
    )
    return kernels

def displacement_scaling(params):
    return params['shear_modulus'] / params['length_scale']

def traction_scaling(params):
    return 1.0

field_types = dict()
field_types[('continuous', 'displacement')] = lambda p: displacement_scaling(p)
field_types[('continuous', 'traction')] = lambda p: traction_scaling(p)
field_types[('discontinuous', 'slip')] = lambda p: displacement_scaling(p)
field_types[('discontinuous', 'crack_traction')] = lambda p: traction_scaling(p)

def get_BIEs(params):
    bies = []
    for name, field_name in field_types.keys():
        if field_name == 'displacement' or field_name == 'slip':
            bies.append(get_displacement_BIE(name, params))
        elif field_name == 'traction' or field_name == 'crack_traction':
            bies.append(get_traction_BIE(name, params))
    return bies

def get_displacement_BIE(obs_mesh_name, params):
    return dict(
        obs_mesh = obs_mesh_name,
        mass_term = dict(
            src_mesh = obs_mesh_name,
            function = 'displacement',
            multiplier = 1.0
        ),
        terms = displacement_BIE_terms(obs_mesh_name, params['gravity']),
        scaling = displacement_scaling(params)
    )

def get_traction_BIE(obs_mesh_name, params):
    return dict(
        obs_mesh = obs_mesh_name,
        mass_term = dict(
            src_mesh = obs_mesh_name,
            function = 'traction',
            multiplier = 1.0
        ),
        terms = traction_BIE_terms(obs_mesh_name, params['gravity']),
        scaling = traction_scaling(params)
    )

def displacement_BIE_terms(obs_mesh_name, gravity):
    terms = [
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'traction',
            function = 'displacement',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'displacement',
            function = 'traction',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'discontinuous',
            kernel = 'traction',
            function = 'slip',
            multiplier = -1.0
        )
    ]
    if gravity:
        terms.append(dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'gravity_displacement',
            function = 'gravity',
            multiplier = 1.0
        ))
    return terms

def traction_BIE_terms(obs_mesh_name, gravity):
    terms = [
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'hypersingular',
            function = 'displacement',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'adjoint_traction',
            function = 'traction',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'discontinuous',
            kernel = 'hypersingular',
            function = 'slip',
            multiplier = -1.0
        )
    ]
    if gravity:
        terms.append(dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'gravity_traction',
            function = 'gravity',
            multiplier = 1.0
        ))
    return terms

