def default_params():
    #TODO: Load this from a file? Or move it to spec?
    return dict(
        obs_order = 3,
        singular_steps = 8,
        far_threshold = 3.0,
        sinh_order = 12,
        solver_tol = 1e-5,
        poisson_ratio = 0.25,
        shear_modulus = 30e9,
        length_scale = 1.0,
        fmm_order = 30, #TODO: fmm_order should be different for each kernel
        dense = False,
        gravity = False
    )

def add_default_parameters(input_params):
    params = default_params()
    for k, v in input_params.iteritems():
        params[k] = v
    return params
