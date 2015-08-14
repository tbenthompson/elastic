import numpy as np

"""
    Defines various information concerning the names of meshes, the fields defined on 
    those meshes and the boundary integral equations that are evaluated on the meshes.
    The "continuous" mesh has displacement continuity from one side to the other while
    the "discontinuous" mesh has a displacement jump from one side to the other.
    "traction" and "displacement" are defined on the continuous mesh and 
    "crack_traction" and "slip" are defined on the discontinuous mesh.
    
    Boundary integral equations are defined as the sum of a bunch of terms
    which look like:
    \int_{S_{obs}} \int_{S_{src}} K(x, y) f(y) dy dx
    The term specifications provide a obs_mesh, a src_mesh, a kernel function, 
    the function f(y), and a constant multiplier.
"""

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

def residual_bie_from_field_name(mesh_name, field_name, params):
    if field_name == 'displacement':
        return get_displacement_BIE(mesh_name, 'displacement', 'displacement', params)
    elif field_name == 'traction':
        return get_traction_BIE(mesh_name, 'traction', 'traction', params)
    elif field_name == 'crack_traction':
        return get_displacement_BIE(mesh_name, 'slip', 'crack_traction', params)
    elif field_name == 'slip':
        return get_displacement_BIE(mesh_name, 'slip', 'slip', params)
    else:
        return 'not a valid bie'

def get_residual_BIEs(params):
    bies = []
    for mesh_name, field_name in field_types:
        bies.append(residual_bie_from_field_name(mesh_name, field_name, params))
    return bies

def get_displacement_BIE(obs_mesh_name, displacement_field, unknown_field, params):
    terms = displacement_BIE_terms(obs_mesh_name, params['gravity'])
    for t in terms:
        t['unknown_field'] = unknown_field
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = unknown_field,
        mass_term = dict(
            src_mesh = obs_mesh_name,
            function = displacement_field,
            multiplier = -1.0
        ),
        terms = terms
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
    terms = traction_BIE_terms(obs_mesh_name, params['gravity'])
    for t in terms:
        t['unknown_field'] = unknown_field
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = unknown_field,
        mass_term = dict(
            src_mesh = obs_mesh_name,
            function = traction_field,
            multiplier = 1.0
        ),
        terms = terms
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

