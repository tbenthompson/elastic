import numpy as np
from elastic import Element, execute
from elastic.solver import form_system, dense_matrix, dense_rhs
import matplotlib.pyplot as plt

def solve_and_plot(es, params, L):
    result = execute(2, es, params)
    x, y = np.meshgrid(np.linspace(0, L, 20), np.linspace(0, L, 20))
    pts = np.array([x.flatten(), y.flatten()]).T
    disp = result.interior_displacement(pts)
    plt.quiver(pts[:, 0], pts[:, 1], disp[0], disp[1])
    plt.xlim([-0.1 * L, 1.1 * L])
    plt.ylim([-0.1 * L, 1.1 * L])
    plt.show()

    f = result.input.meshes['traction'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    s = result.soln[('traction', 'displacement')]
    plt.plot(xs, s[0])
    plt.plot(xs, s[1])
    plt.show()

def box_mesh(L, origin = (0, 0)):
    refine = 4
    es = []
    es.append(Element(
        [[origin[0], origin[1] + L], [origin[0], origin[1]]],
        [[0, 0], [0, 0]], "displacement", refine
    ))
    es.append(Element(
        [[origin[0], origin[1]], [origin[0] + L, origin[1]]],
        [[0, 0], [0, 0]], "displacement", refine
    ))
    es.append(Element(
        [[origin[0] + L, origin[1]], [origin[0] + L, origin[1] + L]],
        [[0, 0], [0, 0]], "displacement", refine
    ))
    es.append(Element(
        [[origin[0] + L, origin[1] + L], [origin[0], origin[1] + L]],
        [[0, 1], [0, 1]], "traction", refine
    ))
    return es

def build_system(L, obs_order, origin = (0, 0)):
    es = box_mesh(L, origin = origin)
    params = dict(
        shear_modulus = 1.0,
        dense = True,
        singular_steps = 8,
        obs_order = obs_order,
        sinh_order = 12,
        length_scale = L
    )
    # solve_and_plot(es, params, L)
    internal_data = form_system(2, es, params)
    matrix = dense_matrix(*internal_data)
    # print(np.linalg.cond(matrix))
    rhs = dense_rhs(*internal_data)
    return matrix, rhs

def obs_order_refine():
    L = 500
    es = box_mesh(L)
    for i in range(15):
        obs_order = 3 + 10 * (i - 1)
        if i == 0:
            obs_order = 3
        if i == 1:
            obs_order = 4
        matrix, rhs = build_system(L, obs_order)
        if i >= 1:
            l2_error = np.sqrt(np.mean((matrix - old_matrix) ** 2) / matrix.size)
            l2_matrix = np.sqrt(np.mean((matrix) ** 2) / matrix.size)
            linf_error = np.max(np.abs(matrix - old_matrix))
            linf_matrix = np.max(np.abs(matrix))
            print('For order: ' + str(obs_order))
            print('L2 error: ' + str(l2_error / l2_matrix))
            print('Linf error: ' + str(linf_error / linf_matrix))
            error = np.abs(matrix - old_matrix)
            plt.imshow(error, interpolation = 'none')
            plt.colorbar()
            plt.show()
        old_matrix = matrix

def size_change():
    matrix1, rhs1 = build_system(1, 5, (0, 0))
    matrix3, rhs3 = build_system(100, 5, (0, 0))
    matrix4, rhs4 = build_system(100, 50, (0, 0))
    import ipdb; ipdb.set_trace()

def helper(L):
    es = box_mesh(L)
    params = dict(
        shear_modulus = 1.0,
        dense = True,
        singular_steps = 8,
        obs_order = 3,
        sinh_order = 50,
        length_scale = L
    )

    #TODO: look at the Linf convergence of the matrix.
    #TODO: can i make a systematic tools to do this?

    result = execute(2, es, params)

    f = result.input.meshes['displacement'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    s = result.soln[('displacement', 'traction')]
    plt.quiver(xs / L, ys / L, s[0] / L, s[1] / L)
    plt.show()

    f = result.input.meshes['traction'].facets
    xs = f[:,:,0].reshape(f.shape[0] * f.shape[1])
    ys = f[:,:,1].reshape(f.shape[0] * f.shape[1])
    s = result.soln[('traction', 'displacement')]
    plt.plot(xs / L, s[0] / L)
    plt.plot(xs / L, s[1] / L)
    return xs, s

def size_change2():
    xs, s1 = helper(1)
    for L2log in range(1, 10, 1):
        L2 = 10 ** L2log
        xs, s2 = helper(L2)
        plt.show()
        errx = s1[0] - (s2[0] / L2)
        erry = s1[1] - (s2[1] / L2)
        print(L2)
        print(np.max(np.abs(errx)))
        print(np.max(np.abs(erry)))
        # matrix, rhs = build_system(L2, 3, (0, 0))
        # print(rhs)
        # plt.plot(xs / L2, errx)
        # plt.plot(xs / L2, erry)
        # plt.show()

if __name__ == '__main__':
    size_change2()

def jagged_mesh():
    refine = 4
    es = []
    es.append(Element([[0, 1], [0, 0]], [[0, 0], [0, 0]], "displacement", refine))
    es.append(Element([[0, 0], [0.25, -1]], [[-0.1, 0], [-0.1, 0]], "traction", refine))
    es.append(Element([[0.25, -1], [0.5, 0]], [[0.0, 0], [0.0, 0]], "displacement", refine))
    es.append(Element([[0.5, 0], [0.75, -1]], [[-0.1, 0], [-0.1, 0]], "traction", refine))
    es.append(Element([[0.75, -1], [1, 0]], [[-0.1, 0], [-0.1, 0]], "traction", refine))
    es.append(Element([[1, 0], [1, 1]], [[0, 0], [0, 0]], "displacement", refine))
    es.append(Element([[1, 1], [0.75, 2]], [[0.1, 0], [0.1, 0]], "traction", refine))
    es.append(Element([[0.75, 2], [0.5, 1]], [[0.0, 0], [0.0, 0]], "displacement", refine))
    es.append(Element([[0.5, 1], [0.25, 2]], [[0.1, 0], [0.1, 0]], "traction", refine))
    es.append(Element([[0.25, 2], [0, 1]], [[0.1, 0], [0.1, 0]], "traction", refine))
    return es
