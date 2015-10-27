from elastic import *
from elastic.interior_mesh_builder import build_interior_mesh

from errors import check_error, check_interior_error
from stress_fnc import *
import numpy as np
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)


def build_stress_fnc(rho, g, k1, k2, k3, k4):
    def stress_fnc(x, y):
        return rho * g / 2.0 * (k1 * x ** 2 * y + k2 * x * y ** 2) + \
               rho * g / 6.0 * (k3 * x ** 3 + k4 * y ** 3)
    return stress_fnc

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

def basal_friction_builder(calc_traction, mu_b):
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
    alpha = np.radians(0.0) # surface slope
    beta = np.radians(5.0) # angle of basal thrust
    theta = alpha + beta # taper

    x0 = 100 * 1000.0 # width of wedge (m)
    def basal_depth(x):
        return x * np.tan(theta)
    d0 = basal_depth(x0) # depth at backstop (m)

    rho = 1500 #density
    g = 9.8 #strength of gravity
    gravity_vector = np.array([
        -rho * g * np.sin(alpha),
        rho * g * np.cos(alpha)
    ])

    mu_b = 0.25 # basal friction
    a = 0.05

    k1 = 0.0
    k2 = (-a * (1 - mu_b * np.tan(theta)) + np.cos(alpha)) / np.tan(theta)
    k3 = 0.0
    k4 = ((a * (3 * (1 - mu_b * np.tan(theta)) * np.cos(theta) ** 2 - 1))
            / (np.sin(theta) ** 2)) - \
        (2 * np.cos(alpha) / (np.tan(theta) ** 2)) - \
        (np.sin(alpha) / np.tan(theta))
    print((k1, k2, k3, k4))

    stress_fnc = build_stress_fnc(rho, g, k1, k2, k3, k4)
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
    fric_bc = basal_friction_builder(trac_calc, mu_b)

    refine = 0
    pts = [[0, 0], [x0, 0], [x0, d0]]
    pts.append(pts[0])
    es = []
    es.extend(line([pts[0], pts[1]], refine, trac_bc))
    es.extend(line([pts[1], pts[2]], refine, trac_bc))
    es.extend(line([pts[2], pts[3]], refine, fric_bc))

    # problems are partially remedied by increasing the src_far_order or refining
    # both solver and interior evaluation are affected
    # check the funny points along the right edge
    # using sinh quadrature everywhere seems to fix a lot of the problems
    # using adaptive quadrature makes everything go haywire!
    # theory:
    # quadrature problems derive from the failure to precisely cancel singularities
    # at the junction between elements
    # theory2:
    # something about the choice of quadrature method depending on the distance from
    # an element. In particular, the use of limits?
    # I need to do a test with a very simple triangle... Three sides, one much shorter than the others... that's what I'm doing, but how can I simplify it? Focus in on one element?

    # zoom in on 98000-100000, 0-800
    params = dict(
        dense = True,
        length_scale = 1.0,
        obs_near_order = 5,
        obs_far_order = 5,
        src_far_order = 12,
        far_threshold = 1000000000000.0,
        gravity = True,
        gravity_vector = gravity_vector
    )

    import elastic.interface

    executor = elastic.interface.Executor(2, es, params)
    system = executor.assemble(elastic.bie_spec.get_BIEs(executor.params))
    dense_solver = elastic.dense_solver.DenseSolver(executor.tbem, executor.params)
    A = dense_solver.get_A(executor.dof_map, executor.constraint_matrix, system)
    b = dense_solver.get_b(executor.dof_map, executor.constraint_matrix, system)
    sing_vecs = np.linalg.svd(A)[1]
    import ipdb; ipdb.set_trace()

    # filename = 'huwangsave'
    result = execute(2, es, params)
    # result.save(filename)
    # result = Result.load(filename)


    disp = result.soln[('continuous', 'displacement')]
    trac = result.soln[('continuous', 'traction')]
    fs = result.meshes['all_mesh'].facets
    xs = fs.reshape((fs.shape[0] * fs.shape[1], fs.shape[2]))
    import matplotlib.pyplot as plt
    plt.plot(disp[0], 'b-.')
    plt.plot(disp[1], 'b-')
    plt.show()
    plt.plot(trac[0], 'k-.')
    plt.plot(trac[1], 'k-')
    plt.show()
    plt.quiver(xs[:, 0], xs[:, 1], disp[0], disp[1])
    plt.show()
    plt.quiver(xs[:, 0], xs[:, 1], trac[0], trac[1])
    plt.show()

    check_error(result, 'continuous', 'traction', trac_calc, 1e-4)

    mesh = build_interior_mesh(
        result.meshes['all_mesh'].facets
    ).refine(2e6)

    interior_pts = mesh.pts

    # interior_pts = np.array([[pts[1][0] * 0.999, pts[2][1] * 0.999]])
    normals_x = np.array([
        np.ones(interior_pts.shape[0]), np.zeros(interior_pts.shape[0])
    ]).T
    normals_y = np.array([
        np.zeros(interior_pts.shape[0]), np.ones(interior_pts.shape[0])
    ]).T
    u_interior = result.interior_displacement(interior_pts)
    tx_interior = result.interior_traction(interior_pts, normals_x)
    ty_interior = result.interior_traction(interior_pts, normals_y)
    fields = [
        u_interior[0], u_interior[1],
        tx_interior[0], tx_interior[1],
        ty_interior[0], ty_interior[1]
    ]
    names = ['ux', 'uy', 'sxx', 'sxy', 'syx', 'syy']

    import matplotlib.pyplot as plt
    for component,name in zip(fields, names):
        plt.figure()
        plt.tricontourf(
            interior_pts[:, 0], interior_pts[:, 1], mesh.tris, component,
            levels = np.linspace(np.min(component), np.max(component), 30),
            extend = 'both'
        )
        plt.colorbar()
        plt.xlabel('$x$ (m)')
        plt.ylabel('$y$ (m)')
        plt.title(name)
    plt.show()

    check_interior_error(interior_pts, normals_x, fields[2:], calc_stress, 3e-3)

if __name__ == '__main__':
    test_hu_wang()
