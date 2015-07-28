from elastic import *

from errors import check_error, check_interior_error
from stress_fnc import *
import numpy as np
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

alpha = np.radians(3.5) # surface slope
beta = np.radians(10.0) # angle of basal thrust
theta = alpha + beta # taper

x0 = 100 * 1000.0 # width of wedge (m)

def basal_depth(x):
    return x * np.tan(theta)
d0 = basal_depth(x0) # depth at backstop (m)

rho = 2500 #density
g = 9.8 #strength of gravity
gravity_vector = np.array([
    -rho * g * np.sin(alpha),
    rho * g * np.cos(alpha)
])

mu_b = 0.7 # basal friction

a = 0.5
k1 = 0.0
k2 = (-a * (1 - mu_b * np.tan(theta)) + np.cos(alpha)) / np.tan(theta)
k3 = 0.0
k4 = ((a * (3 * (1 - mu_b * np.tan(theta)) * np.cos(theta) ** 2 - 1))
        / (np.sin(theta) ** 2)) - \
    (2 * np.cos(alpha) / (np.tan(theta) ** 2)) - \
    (np.sin(alpha) / np.tan(theta))
print((k1, k2, k3, k4))

def stress_fnc(x, y):
    return rho * g / 2.0 * (k1 * x ** 2 * y + k2 * x * y ** 2) + \
           rho * g / 6.0 * (k3 * x ** 3 + k4 * y ** 3)

def basal_friction(pts, normal_bc, friction_coefficient):
    from elastic.constraints import Constraint, BCConstraint, Term
    dim = len(pts)
    if dim == 3:
        return 'unimplemented!'
    cs = [
        Constraint([
            Term('traction', 'normal', friction_coefficient),
            Term('traction', 'tangential0', -1.0),
        ], [0, 0]),
        BCConstraint('traction', 'normal', normal_bc)
    ]
    return dict(pts = pts, type = 'continuous', constraints = cs)

def basal_friction_builder(calc_traction):
    def bc(pts):
        nx = -(pts[1, 1] - pts[0, 1])
        ny = pts[1, 0] - pts[0, 0]
        nmag = np.sqrt(nx ** 2 + ny ** 2)
        nx /= nmag
        ny /= nmag
        t0 = calc_traction(pts[0, :], [nx, ny])
        tn0 = t0[0] * nx + t0[1] * ny
        t1 = calc_traction(pts[1, :], [nx, ny])
        tn1 = t1[0] * nx + t1[1] * ny
        return basal_friction(pts, [tn0, tn1], -mu_b)
    return bc

def test_hu_wang():
    stress_lambdas = form_stress_lambdas(stress_fnc, gravity_vector)
    calc_stress = calc_stress_builder(stress_lambdas)
    trac_calc = calc_traction_builder(calc_stress)

    npts = 200
    x = np.linspace(0, x0, npts)
    y = basal_depth(x)
    nx = -(y[1] - y[0])
    ny = x[1] - x[0]
    nmag = np.sqrt(nx ** 2 + ny ** 2)
    nx /= nmag
    ny /= nmag
    tracs = [trac_calc([x[i], y[i]], [nx, ny]) for i in range(npts)]
    tn = [t[0] * nx + t[1] * ny for t in tracs]
    ts = [-t[0] * ny + t[1] * nx for t in tracs]
    np.testing.assert_almost_equal(ts[1] / tn[1], mu_b, 6)


    trac_bc = traction_bc_builder(trac_calc)
    fric_bc = basal_friction_builder(trac_calc)

    refine = 5
    pts = [[0, 0], [x0, 0], [x0, d0], [0, 0]]
    es = []
    es.extend(line([pts[0], pts[1]], refine, trac_bc))
    es.extend(line([pts[1], pts[2]], refine, trac_bc))
    es.extend(line([pts[2], pts[3]], refine, fric_bc))
    # es.extend(line([pts[2], pts[3]], refine, trac_bc))

    params = dict(
        dense = True,
        singular_steps = 6,
        obs_order = 10,
        sinh_order = 10,
        gravity = True,
        gravity_vector = gravity_vector
    )
    result = execute(2, es, params)
    result.params['obs_order'] = 15

    check_error(result, 'continuous', 'traction', trac_calc, 1e-2)

    n_interior = 80
    x_hat, y_hat = np.meshgrid(
        np.linspace(0, 1, n_interior), np.linspace(0, 1, n_interior)
    )

    x_scale = pts[1][0]
    y_scale = pts[2][1]
    x = x_hat * x_scale
    y = y_hat * (x / x_scale) * y_scale
    pts = np.array([x.flatten(), y.flatten()]).T
    normals_x = np.array([np.ones(pts.shape[0]), np.zeros(pts.shape[0])]).T
    normals_y = np.array([np.zeros(pts.shape[0]), np.ones(pts.shape[0])]).T
    u_interior = result.interior_displacement(pts)
    tx_interior = result.interior_traction(pts, normals_x)
    ty_interior = result.interior_traction(pts, normals_y)
    fields = [
        u_interior[0], u_interior[1],
        tx_interior[0], tx_interior[1],
        ty_interior[0], ty_interior[1]
    ]

    x_mat = pts[:, 0].reshape((n_interior, n_interior))
    y_mat = pts[:, 1].reshape((n_interior, n_interior))
    import matplotlib.pyplot as plt
    for s_component in fields:
        s_mat = s_component.reshape(x_mat.shape)
        plt.figure()
        plt.contourf(
            x_mat, y_mat, s_mat,
            levels = np.linspace(np.min(s_mat), np.max(s_mat), 30)
        )
        plt.colorbar()
    plt.show()

    check_interior_error(pts, normals_x, fields[2:], calc_stress, 1e-3)
