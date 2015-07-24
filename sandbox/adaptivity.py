from elastic import adaptive_execute, line, displacement, traction, mixed, free_slip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import logging
import copy

def adaptive():
    logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

    es = line([[1, 1], [0.5, 1]], 0,
        lambda pts: mixed(pts, ['traction', 'displacement'], [[0, -1], [0, -1]])
    )
    es += line([[0.5, 1], [0.0, 1]], 0,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    )
    es += line([[0.0, 1], [0.0, 0.0]], 0,
        lambda pts: displacement(pts, [[1, 0], [1, 0]])
    )
    es += line([[0.0, 0], [1.0, 0]], 0,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    )
    es += line([[1.0, 0], [1.0, 1]], 0,
        lambda pts: traction(pts, [[0, 0], [0, 0]])
    )
    # es += line([[0.5, 0.5], [0.0, 1]], 0,
    #     lambda pts: free_slip(pts, [0, 0], [[0, 0]])
    # )
    params = dict()
    params['error_threshold'] = 0.005
    params['max_iters'] = 25
    params['refine_fraction'] = 0.2
    params['dense'] = True
    result = adaptive_execute(2, es, params, callback = callback)
    animate()

storage = []

def callback(executor, result):
    # make_plot(executor, result)
    # plt.show()
    storage.append([copy.copy(executor), copy.copy(result)])

def make_plot(executor, result):
    meshes = executor.meshes
    f = meshes['all_mesh'].facets
    plt.clf()
    for i in range(f.shape[0]):
        plt.plot(f[i, :, 0], f[i, :, 1], 'k-', alpha = 0.3)

    ux = result.soln[('continuous', 'displacement')][0]
    uy = result.soln[('continuous', 'displacement')][1]
    squish_factor = 0.25

    f = result.meshes['continuous'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    x_def = xs + squish_factor * ux
    y_def = ys + squish_factor * uy
    plt.plot(x_def, y_def, 'kH-', linewidth = 0.5, markersize = 3)

    f = result.meshes['discontinuous'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    interior_pts = np.array([xs.flatten(), ys.flatten()]).T
    ux, uy = result.interior_displacement(interior_pts)
    x_def = xs + squish_factor * ux
    y_def = ys + squish_factor * uy
    plt.plot(x_def, y_def, 'kH-', linewidth = 0.5, markersize = 3)


    nx = 100
    ny = 100
    x_vec = np.linspace(0, 1, nx)
    y_vec = np.linspace(0, 1, ny)
    x_mat, y_mat = np.meshgrid(x_vec, y_vec)
    interior_pts = np.array([x_mat.flatten(), y_mat.flatten()]).T
    ux, uy = result.interior_displacement(interior_pts)
    ux = ux.reshape(x_mat.shape)
    uy = uy.reshape(y_mat.shape)
    x_mat_def = x_mat + squish_factor * ux
    y_mat_def = y_mat + squish_factor * uy
    plt.contourf(
        x_mat_def, y_mat_def, np.sqrt(ux ** 2 + uy ** 2)
    )
    plt.colorbar()
    plt.axis([-0.1, 1.5, -0.5, 1.3])

    # normals_x = np.array(
    #     [np.ones(interior_pts.shape[0]), np.zeros(interior_pts.shape[0])]
    # ).T
    # sxx, sxy = result.interior_traction(interior_pts, normals_x)
    # sxy = sxy.reshape(x_mat.shape)

def animate():
    def f(step_idx):
        make_plot(*storage[step_idx])
    fig = plt.figure(figsize = (10,10))
    anim = animation.FuncAnimation(fig, f, frames = len(storage))
    anim.save('elastic_refine.gif', writer = 'imagemagick', fps = 1)

if __name__ == '__main__':
    adaptive()
