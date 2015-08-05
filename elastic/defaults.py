import time

def default_params():
    #TODO: Load this from a file? Or move it to spec?
    return dict(
        poisson_ratio = 0.25,
        shear_modulus = 30e9,
        length_scale = 1.0,

        obs_near_order = 4,
        obs_far_order = 4,
        src_far_order = 3,
        far_threshold = 3.0,

        singular_steps = 8,
        sinh_order = 12,
        fmm_order = 30, #TODO: fmm_order should be different for each kernel

        solver_tol = 1e-5,
        dense = False,
        check_condition_number = True,

        gravity = False,
        gravity_vector = [0.0, -9.8 * 2700],

        timing = False,
        timer = time,
    )

def add_default_parameters(input_params):
    params = default_params()
    params.update(input_params)
    return params
