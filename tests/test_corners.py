import elastic

import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

def corners():
    refine = 5

    es = []

    #lower edge is displacement upwards
    es.extend(elastic.line([[0, 0], [1, 0]], refine,
        lambda pts: elastic.displacement(pts, [[1, 0], [1, 0]])
    ))

    #left and right edges are traction free
    es.extend(elastic.line([[1, 0], [1, 2]], refine,
        lambda pts: elastic.traction(pts, [[0, 0], [0, 0]])
    ))
    es.extend(elastic.line([[0, 2], [0, 0]], refine,
        lambda pts: elastic.traction(pts, [[0, 0], [0, 0]])
    ))

    es.extend(elastic.line([[1, 2], [0, 2]], refine,
        lambda pts: elastic.displacement(pts, [[0, 0], [0, 0]])
    ))

    def free_slip_bc(pts):
        return dict(pts = pts, type = 'discontinuous', constraints = [
            elastic.constraints.BCConstraint('crack_traction', 'tangential0', [0, 0]),
            elastic.constraints.BCConstraint('slip', 'normal', [0.0, 0.0])
        ])

    #upper edge is serrated and traction free
    n_corners = 32
    corner_height = 1.01
    corner_refine = int(refine - np.log2(n_corners))
    xs = np.linspace(0, 1, 2 * n_corners + 1)
    for i in range(n_corners):
        pt0 = [xs[2 * i], 1]
        pt1 = [xs[2 * i + 1], corner_height]
        pt2 = [xs[2 * i + 2], 1]
        es.extend(elastic.line([pt1, pt0], corner_refine, free_slip_bc))
        es.extend(elastic.line([pt2, pt1], corner_refine, free_slip_bc))

    params = dict(
        dense = True,
        shear_modulus = 1.0,
        poisson_ratio = 0.25,
        src_far_order = 10,
    )

    result = elastic.execute(2, es, params);

    f = result.meshes['continuous'].facets
    bdry_pts = f.reshape((f.shape[0] * f.shape[1], f.shape[2]))
    which = np.logical_and(bdry_pts[:, 0] < 0.02, bdry_pts[:, 1] > 0.9)
    t = result.soln[('continuous', 'traction')]
    print(np.max(t))
    d = result.soln[('continuous', 'displacement')]
    plt.figure()
    plt.quiver(bdry_pts[:, 0], bdry_pts[:, 1], t[0][:], t[1][:])
    plt.figure()
    plt.quiver(bdry_pts[:, 0], bdry_pts[:, 1], d[0][:], d[1][:])
    plt.show()

    bdry_facets = result.meshes['all_mesh'].facets
    # corners = elastic.interior_mesh_builder.extents_to_box_2d(
    #     [0, 0.99], [0.1, 1.01]
    # )
    # bdry_facets = np.vstack((bdry_mesh, corners))
    mesh = elastic.interior_mesh_builder.build_interior_mesh(bdry_facets)
    mesh.region_plot()
    mesh = mesh.get_regions([mesh.get_region_id_from_pt([0.01, 1.0 - 1e-8])])
    max_tri_area = 0.00002
    mesh = mesh.refine(max_tri_area)
    # mesh.region_plot()

    s_interior = np.array(result.interior_stress(mesh.pts))

    field = np.log(np.abs(s_interior[0, 1, :] + 1e-8))
    plt.figure()
    plt.tripcolor(
        mesh.pts[:, 0], mesh.pts[:, 1], mesh.tris, field, shading = 'gouraud',
        # vmin = -5.0, vmax = 5.0
    )
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    corners()
