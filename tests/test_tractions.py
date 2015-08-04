from elastic import execute, line, traction
from errors import check_error, check_interior_error
import numpy as np
from stress_fnc import *

#stress function =
#a / 6 * x^3 + b / 2 * x^2 * y + c / 2 * x * y^2 + d / 6 * y^3

# TODO: there are some blips and uglies when evaluating the interior tractions
# these disappear when using high gaussian quad orders or refining the mesh
# sufficiently.
# Seems like 4 src quadrature points is sufficient. I think I should decouple
# the source quadrature points from the observation quadrature points again.
# For some reason, there isn't a very good response to refinement which makes
# me think that the problem is in the quadrature.
# The sinh quadrature seems to do just fine!
# The limit quadrature also seems to be just fine!
# The fact that the quadrature allows something like 1% error is good enough for
# the moment, but should be fixed eventually.
# TODO: Other interesting observation. There seems to be a hard lower bound of 1.6s
# on the running time of this file. Weird. Where's it going?

def stress_fnc(x, y):
    a,b,c,d = (0.2, 0.8, -0.7, 0.3)#np.random.rand(4)
    return a / 6 * x ** 3 + b / 2 * x ** 2 * y + c / 2 * x * y ** 2 + d / 6 * y ** 3


def test_tractions():
    calc_stress = calc_stress_builder(form_stress_lambdas(stress_fnc, [0, 0]))
    trac_calc = calc_traction_builder(calc_stress)
    trac_bc = traction_bc_builder(trac_calc)

    refine = 7
    pts = [[0, 0], [100, -100], [100, 0], [0, 0]]
    es = []
    for i in range(len(pts) - 1):
        es.extend(line([pts[i], pts[i + 1]], refine, trac_bc))

    params = dict(
        dense = True,
        shear_modulus = 1.0,
        poisson_ratio = 0.3,
        singular_steps = 8,
        obs_far_order = 12,
        obs_near_order = 12,
        src_far_order = 7,
        sinh_order = 10
    )
    result = execute(2, es, params)

    check_error(result, 'continuous', 'traction', trac_calc, 1e-2)

    n_interior = 40
    x_hat, y_hat = np.meshgrid(
        np.linspace(0, 1, n_interior), np.linspace(0, 1, n_interior)
    )

    x_scale = pts[1][0]
    y_scale = pts[1][1]
    x = x_hat * x_scale
    y = y_hat * (x / x_scale) * y_scale
    pts = np.array([x.flatten(), y.flatten()]).T
    normals_x = np.array([np.ones(pts.shape[0]), np.zeros(pts.shape[0])]).T
    normals_y = np.array([np.zeros(pts.shape[0]), np.ones(pts.shape[0])]).T
    tx_interior = result.interior_traction(pts, normals_x)
    ty_interior = result.interior_traction(pts, normals_y)
    stress_interior = [tx_interior[0], tx_interior[1], ty_interior[0], ty_interior[1]]
    check_interior_error(pts, normals_x, stress_interior, calc_stress, 1e-3)

