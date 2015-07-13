import numpy as np

mesh_types = [
    'continuous',
    'discontinuous'
]

field_types = [
    ('continuous', 'traction'),
    ('continuous', 'displacement'),
    ('discontinuous', 'crack_traction'),
    ('discontinuous', 'slip'),
]

unknowns_to_knowns = dict(
    traction = 'displacement',
    displacement = 'traction',
    crack_traction = 'slip',
    slip = 'crack_traction'
)

field_units = dict(
    displacement = 'displacement',
    slip = 'displacement',
    traction = 'traction',
    crack_traction = 'traction'
)

dimensionless_scaling = dict()
dimensionless_scaling['displacement'] = lambda p: p['shear_modulus']
dimensionless_scaling['traction'] = lambda p: p['length_scale']

solution_scaling = dict()
for k in field_types:
    solution_scaling[k] = dimensionless_scaling[field_units[k[1]]]

integral_scaling = dict()
for k in field_types:
    integral_scaling[k] = dimensionless_scaling[unknowns_to_knowns[field_units[k[1]]]]

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

def bie_from_field_name(mesh_name, field_name, params):
    if field_name == 'displacement':
        return get_traction_BIE(mesh_name, 'traction', 'displacement', params)
    elif field_name == 'traction':
        return get_displacement_BIE(mesh_name, 'displacement', 'traction', params)
    elif field_name == 'crack_traction':
        return get_traction_BIE(mesh_name, 'crack_traction', 'crack_traction', params)
    elif field_name == 'slip':
        return get_traction_BIE(mesh_name, 'crack_traction', 'slip', params)
    else:
        return 'not a valid bie'

def get_BIEs(params):
    bies = []
    for mesh_name, field_name in field_types:
        bies.append(bie_from_field_name(mesh_name, field_name, params))
    return bies

def get_displacement_BIE(obs_mesh_name, displacement_field, unknown_field, params):
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = unknown_field,
        mass_term = dict(
            src_mesh = obs_mesh_name,
            function = displacement_field,
            multiplier = -1.0
        ),
        terms = displacement_BIE_terms(obs_mesh_name, params['gravity'])
    )

def displacement_BIE_terms(obs_mesh_name, gravity):
    terms = [
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'traction',
            function = 'displacement',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'displacement',
            function = 'traction',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'discontinuous',
            kernel = 'traction',
            function = 'slip',
            multiplier = 1.0
        )
    ]
    if gravity:
        terms.append(dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'continuous',
            kernel = 'gravity_displacement',
            function = 'ones',
            multiplier = -1.0
        ))
    return terms

def get_traction_BIE(obs_mesh_name, traction_field, unknown_field, params):
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = unknown_field,
        mass_term = dict(
            src_mesh = obs_mesh_name,
            function = traction_field,
            multiplier = 1.0
        ),
        terms = traction_BIE_terms(obs_mesh_name, params['gravity'])
    )

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
            function = 'ones',
            multiplier = 1.0
        ))
    return terms

