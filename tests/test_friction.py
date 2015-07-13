from elastic import displacement, traction, static_friction, line, execute,\
    free_slip, mixed, slip
import numpy as np

def test_friction():
    refine = 6

    es = []
    es.extend(line([[1, 1], [0, 1]], refine,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    ))
    es.extend(line([[0, 1], [0, 0]], refine,
        lambda pts: mixed(pts, ['displacement', 'traction'], [[1, 0], [1, 0]])
    ))
    es.extend(line([[0, 0], [1, 0]], refine,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    ))
    es.extend(line([[1, 0], [1, 1]], refine,
        lambda pts: mixed(pts, ['displacement', 'traction'], [[0, 0], [0, 0]])
    ))
    es.extend(line([[0.2, 0.3], [0.553182, 1.0]], refine,
        lambda pts: free_slip(pts, [0, 0], [[0, 0]])
    ))

    result = execute(2, es, dict(dense = True))

    plotter(result, 'continuous', 'displacement', 0)

    nx = 100
    ny = 100
    x_vec = np.linspace(0, 1, nx)
    y_vec = np.linspace(0, 1, ny)
    x_mat, y_mat = np.meshgrid(x_vec, y_vec)
    interior_pts = np.array([x_mat.flatten(), y_mat.flatten()]).T
    normals_x = np.array(
        [np.ones(interior_pts.shape[0]), np.zeros(interior_pts.shape[0])]
    ).T
    normals_y = np.array(
        [np.zeros(interior_pts.shape[0]), np.ones(interior_pts.shape[0])]
    ).T
    ux, uy = result.interior_displacement(interior_pts)
    sxx, sxy = result.interior_traction(interior_pts, normals_x)
    sxy, syy = result.interior_traction(interior_pts, normals_y)

    ux_mat = ux.reshape(x_mat.shape)
    uy_mat = uy.reshape(x_mat.shape)
    sxx_mat = sxx.reshape(x_mat.shape)
    sxy_mat = sxy.reshape(x_mat.shape)
    syy_mat = syy.reshape(x_mat.shape)

    import matplotlib.pyplot as plt
    plt.contourf(x_mat, y_mat, np.sqrt(ux_mat ** 2 + uy_mat ** 2))
    plt.colorbar()
    def skip_and_make_vec(f, k):
        skipped = f[::k,::k]
        return skipped.reshape(skipped.size)
    plt.quiver(
        skip_and_make_vec(x_mat, 5), skip_and_make_vec(y_mat, 5),
        skip_and_make_vec(ux_mat, 5), skip_and_make_vec(uy_mat, 5)
    )
    plt.show()

def plotter(result, mesh_name, field_name, x_component):
    f = result.meshes[mesh_name].facets
    xs = f[:, :, x_component].reshape(f.shape[0] * f.shape[1])
    data = result.soln[(mesh_name, field_name)][0]
    import matplotlib.pyplot as plt
    plt.plot(xs, data)
    plt.show()

