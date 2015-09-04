import elastic

import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

def test_corners():
    refine = 5

    es = []

    #lower edge is displacement upwards
    es.extend(elastic.line([[0, 0], [1, 0]], refine,
        lambda pts: elastic.displacement(pts, [[0, 1], [0, 1]])
    ))

    #left and right edges are traction free
    es.extend(elastic.line([[1, 0], [1, 1]], refine,
        lambda pts: elastic.traction(pts, [[0, 0], [0, 0]])
    ))
    es.extend(elastic.line([[0, 1], [0, 0]], refine,
        lambda pts: elastic.traction(pts, [[0, 0], [0, 0]])
    ))

    #upper edge is serrated and traction free
    n_corners = 32
    corner_height = 1.01
    corner_refine = int(refine - np.log2(n_corners))
    xs = np.linspace(0, 1, 2 * n_corners + 1)
    for i in range(n_corners):
        pt0 = [xs[2 * i], 1]
        pt1 = [xs[2 * i + 1], corner_height]
        pt2 = [xs[2 * i + 2], 1]
        es.extend(elastic.line([pt0, pt1], corner_refine,
            lambda pts: elastic.traction(pts, [[0, 0], [0, 0]])
        ))
        es.extend(elastic.line([pt1, pt2], corner_refine,
            lambda pts: elastic.traction(pts, [[0, 0], [0, 0]])
        ))

    params = dict(
        dense = True,
        shear_modulus = 1.0,
        poisson_ratio = 0.25
    )

    result = elastic.execute(2, es, params);

    bdry_mesh = result.meshes['all_mesh'].facets
    corners = elastic.interior_mesh_builder.extents_to_box_2d([0, 0.995], [1, 1.02])
    bdry_mesh = np.vstack((bdry_mesh, corners))
    import tbempy
    mesh = elastic.interior_mesh_builder.build_interior_mesh(
        tbempy.TwoD.Mesh(bdry_mesh)
    )
    mesh = mesh.get_region(mesh.get_region_id_from_pt([0.5, 1.0 - 1e-8]))
    max_tri_area = 0.0000002
    mesh = mesh.refine(max_tri_area)
    mesh.region_plot()

    s_interior = np.array(result.interior_stress(mesh.pts))

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
