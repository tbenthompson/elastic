import tbempy.TwoD as tbem
from elastic import Element, execute
import matplotlib.pyplot as plt
import numpy as np

def test_gravity():
    L = 1000.0
    H = 1000.0
    g = 9.8
    rho = 2700.0
    syy = rho * g * H
    R = 8
    es = []
    es.append(Element([[L, 0], [-L, 0]], [[0, 0], [0, 0]], "traction", R))
    # es.append(Element([[-L, 0], [-L, -H]], [[0, 0], [0, 0]], "displacement", R))
    # es.append(Element([[-L, -H], [L, -H]], [[0, syy], [0, syy]], "traction", R))
    # es.append(Element([[L, -H], [L, 0]], [[0, 0], [0, 0]], "displacement", R))

    es.append(Element([[L, 0], [-L, 0]], [[0, 0], [0, 0]], "gravity", R))
    # es.append(Element([[-L, 0], [-L, -H]], [[0, 0], [0, 0]], "gravity", R))
    # es.append(Element([[-L, -H], [L, -H]], [[0, 0], [0, 0]], "gravity", R))
    # es.append(Element([[L, -H], [L, 0]], [[0, 0], [0, 0]], "gravity", R))
    params = dict(
        gravity = True,
        gravity_vector = [0.0, -rho * g],
        dense = True
    )
    result = execute(2, es, params)

    ptsy = np.linspace(-1.0, -0.05 * L, 1000)
    pts = np.vstack((np.zeros(ptsy.shape[0]), ptsy)).T
    normalsx = np.array([[1, 0]] * ptsy.shape[0])
    normalsy = np.array([[0, 1]] * ptsy.shape[0])
    interior_tracx = result.interior_traction(pts, normalsx)
    interior_tracy = result.interior_traction(pts, normalsy)
    import ipdb; ipdb.set_trace()
    plt.plot(ptsy, (1.0 / 3.0) * rho * g * ptsy, 'b-.', label = 'sxx correct')
    plt.plot(ptsy, rho * g * ptsy, 'k-.', label = 'syy correct')
    plt.plot(ptsy, interior_tracx[0], 'b-', label = 'sxx')
    plt.plot(ptsy, interior_tracx[1], 'r-', label = 'sxy')
    plt.plot(ptsy, interior_tracy[1], 'k-', label = 'syy')
    plt.legend()
    plt.show()

    f_trac = result.input.meshes['traction'].facets
    xs_trac = f_trac[:,:,0].reshape(f_trac.shape[0] * f_trac.shape[1])
    ys_trac = f_trac[:,:,1].reshape(f_trac.shape[0] * f_trac.shape[1])
    f_disp = result.input.meshes['displacement'].facets
    xs_disp = f_disp[:,:,0].reshape(f_disp.shape[0] * f_disp.shape[1])
    ys_disp = f_disp[:,:,1].reshape(f_disp.shape[0] * f_disp.shape[1])
    s_disp = result.soln[('traction', 'displacement')]
    s_trac = result.soln[('displacement', 'traction')]
    # indices = np.logical_and((-H + (H / 10) < ys_disp), (ys_disp  < (-H / 10)))

    # plt.plot(ys_disp[indices], s_trac[0][indices], 'b-')
    # plt.plot(ys_disp[indices], s_trac[1][indices], 'r-')
    # plt.plot(ys_disp[indices], rho * g * (1.0 / 3.0) * ys_disp[indices], 'k-')
    # plt.show()
    # # plt.quiver(xs, ys, s_disp[0], s_disp[1])
    # plt.quiver(xs_disp, ys_disp, s_trac[0], s_trac[1])
    # plt.xlim([-L - (L / 10), L + (L / 10)])
    # plt.ylim([-H - (H / 10), (H / 10)])
    # plt.show()
    plt.quiver(xs_trac, ys_trac, s_disp[0], s_disp[1])
    plt.xlim([-L - (L / 10), L + (L / 10)])
    plt.ylim([-H - (H / 10), (H / 10)])
    plt.show()


if __name__ == '__main__':
    test_gravity()
