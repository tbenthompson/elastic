def default_params():
    #TODO: Load this from a file? Or move it to spec?
    return dict(
        poisson_ratio = 0.25,
        shear_modulus = 30e9,
        length_scale = 1.0,
        obs_near_order = 4,
        obs_far_order = 4,
        src_far_order = 3,
        singular_steps = 8,
        far_threshold = 3.0,
        sinh_order = 12,
        solver_tol = 1e-5,
        fmm_order = 30, #TODO: fmm_order should be different for each kernel
        dense = False,
        gravity = False,
        gravity_vector = [0.0, -9.8 * 2700],
        timing = False
    )

def add_default_parameters(input_params):
    params = default_params()
    params.update(input_params)
    return params
