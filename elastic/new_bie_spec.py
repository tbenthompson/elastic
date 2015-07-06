def displacement_BIE_terms(gravity = False):
    terms = [
        dict(
            obs_mesh = 'continuous',
            src_mesh = 'continuous',
            kernel = 'traction',
            function = 'displacement',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = 'continuous',
            src_mesh = 'continuous',
            kernel = 'displacement',
            function = 'traction',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = 'continuous',
            src_mesh = 'discontinuous',
            kernel = 'traction',
            function = 'slip',
            multiplier = -1.0
        )
    ]
    if gravity:
        terms.append(dict(
            obs_mesh = 'continuous',
            src_mesh = 'continuous',
            kernel = 'gravity_displacement',
            function = 'gravity',
            multiplier = 1.0
        ))
    return terms

def traction_BIE_terms(obs_mesh_name, params):
    terms = [
        dict(
            obs_mesh = 'continuous',
            src_mesh = 'continuous',
            kernel = 'hypersingular',
            function = 'displacement',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = 'continuous',
            src_mesh = 'continuous',
            kernel = 'adjoint_traction',
            function = 'traction',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = 'continuous',
            src_mesh = 'discontinuous',
            kernel = 'hypersingular',
            function = 'slip',
            multiplier = -1.0
        )
    ]
    if gravity:
        terms.append(dict(
            obs_mesh = 'continuous',
            src_mesh = 'gravity',
            kernel = 'gravity_traction',
            function = 'gravity',
            multiplier = 1.0
        ))
    return terms

