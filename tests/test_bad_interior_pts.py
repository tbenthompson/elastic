import elastic

import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

def test_corners():
    refine = 0
    es = [elastic.displacement(np.array([[0, 0], [1, 0]]), [[0, 1], [0, 1]])]
    params = dict(dense = True, shear_modulus = 1.0, poisson_ratio = 0.25)
    result = elastic.execute(2, es, params);

    pts = []
    for e in es:
        pts.append(0.5 * e['pts'][1, :] + 0.5 * e['pts'][0, :])
    pts = np.array(pts)

    stresses = np.array(result.interior_stress(pts))
    assert(not np.any(np.isnan(stresses)))


def print_failures(pts, stresses):
    fails = []
    for i in range(pts.shape[0]):
        if np.isnan(stresses[0, 0, i]):
            fails.append(pts[i, :])
    fails = np.array(fails)

    np.set_printoptions(precision = 17)
    print(fails)

def interior_plot(mesh, s_interior):
    field = s_interior[0, 0, :]
    plt.figure()
    plt.tripcolor(
        mesh.pts[:, 0], mesh.pts[:, 1], mesh.tris, field, shading = 'gouraud',
        vmin = -5.0, vmax = 5.0
    )
    plt.colorbar()
    plt.show()

if __name__ == '__main__':
    test_corners()
