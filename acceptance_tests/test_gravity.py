import tbempy.TwoD as tbem
from elastic import Element, execute
from elastic.mesh_gen import circle
import matplotlib.pyplot as plt
import numpy as np

def test_gravity():
    r = 1000.0
    es = []
    es.extend(circle([0, 0], r, 6, "displacement", lambda x: np.zeros_like(x), False))
    es.extend(circle([0, 0], r, 6, "gravity", lambda x: np.zeros_like(x), False))

    params = dict(
        obs_order = 5,
        singular_steps = 12,
        near_tol = 1e-6,
        shear_modulus = 30e9,
        poisson_ratio = 0.25,
        gravity = True,
        gravity_vector = [0.0, -9.8 * 2700],
        dense = True
    )
    result = execute(2, es, params)

    pts = []
    x = y = np.linspace(-r, r, 100)
    x_mat, y_mat = np.meshgrid(x, y)
    pts = np.array([x_mat.flatten(), y_mat.flatten()]).T
    normalsx = np.array([[1, 0]] * pts.shape[0])
    normalsy = np.array([[0, 1]] * pts.shape[0])

    ux, uy = result.interior_displacement(pts)
    sxx, sxy = result.interior_traction(pts, normalsx)
    sxy, syy = result.interior_traction(pts, normalsy)

    mask = x_mat ** 2 + y_mat ** 2 > r ** 2
    sxy_mat = sxy.reshape(x_mat.shape)
    sxy_mat[mask] = np.nan
    x_mat[mask] = np.nan
    y_mat[mask] = np.nan
    plt.contourf(x_mat, y_mat, sxy_mat,
        levels = np.linspace(-8e7, 8e7, 25), extend = 'both')
    plt.xlim([-r, r])
    plt.ylim([-r, r])
    plt.show()

    import sys; sys.exit()

    r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    plt.figure()
    plt.plot(r, uy, '.')
    plt.title('uy')
    plt.figure()
    plt.plot(pts[:, 1], sxx, '.')
    plt.title('sxx')
    plt.figure()
    plt.plot(pts[:, 0], sxy, '.')
    plt.title('sxy')
    plt.figure()
    plt.plot(pts[:, 1], syy, '.')
    plt.title('syy')
    plt.show()
    # plt.figure()
    # plt.quiver(pts[:, 0], pts[:, 1], ux, uy)

    # plt.figure()
    # plt.quiver(pts[:, 0], pts[:, 1], sxx, sxy)
    # plt.figure()
    # plt.quiver(pts[:, 0], pts[:, 1], sxy, syy)
    # plt.show()

if __name__ == '__main__':
    test_gravity()
