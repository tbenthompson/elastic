from elastic.interior_mesh_builder import *
import numpy as np

def box_mesh():
    return np.array([
        [[0, 0], [1, 0]],
        [[1, 0], [1, 1]],
        [[1, 1], [0, 1]],
        [[0, 1], [0, 0]]
    ])

def test_determine_extents():
    min_corner, max_corner = determine_extents(box_mesh())
    np.testing.assert_almost_equal(min_corner, [0, 0])
    np.testing.assert_almost_equal(max_corner, [1, 1])

def test_determine_extents_negative_quadrant():
    min_corner, max_corner = determine_extents(box_mesh() - 1)
    np.testing.assert_almost_equal(min_corner, [-1, -1])
    np.testing.assert_almost_equal(max_corner, [0, 0])

def test_expand_extents():
    min_corner, max_corner = expand_extents(*determine_extents(box_mesh()))
    np.testing.assert_almost_equal(min_corner, [-0.5, -0.5])
    np.testing.assert_almost_equal(max_corner, [1.5, 1.5])

def test_extents_to_box():
    result = extents_to_box_2d(*determine_extents(box_mesh()))
    np.testing.assert_almost_equal(result, box_mesh())

def test_interior_meshing():
    facets = add_extent_surface(box_mesh())
    mesh = InteriorMeshBuilder(facets)
    assert(len(mesh.tris) == 14)
    assert(not mesh.are_across_an_input_facet(4, 11))
    assert(mesh.are_across_an_input_facet(11, 9))
    regions = mesh.identify_regions()

    assert(np.max(regions) == 1)
    assert(np.all(regions == [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0]))

def test_adjacent_tris():
    tris = np.array([
        [0, 1, 2], [2, 1, 3], [1, 5, 4], [1, 4, 3]
    ])
    pts_to_tris_map = build_pt_to_tri_map(tris)

    tri_indices = find_adjacent_tris_indices(tris, [1, 3, 4], pts_to_tris_map)
    assert(tri_indices == [1, 2])

    tri_indices = find_adjacent_tris_indices(tris, [0, 1, 2], pts_to_tris_map)
    assert(tri_indices == [1])

def test_build_pt_to_tri_map():
    tris = np.array([
        [0, 1, 2], [2, 1, 3], [1, 5, 4], [1, 4, 3]
    ])
    map = build_pt_to_tri_map(tris)
    assert(map[1] == [0, 1, 2, 3])
    assert(map[2] == [0, 1])

def make_pretty_region_plots():
    facets = np.array([
        [[0, 0], [1, 0]],
        [[1, 0], [1, 1]],
        [[1, 1], [0, 1]],
        [[0, 1], [0, 0]],
        [[0, 0], [1, 1]]
    ])
    facets = add_extent_surface(facets)
    facets = np.vstack((facets, [[[-0.5, -0.5], [0, 0]], [[1.0, 1.0], [1.5, 1.5]]]))
    mesh = InteriorMesh(facets)
    mesh.region_plot()

if __name__ == '__main__':
    make_pretty_region_plots()
