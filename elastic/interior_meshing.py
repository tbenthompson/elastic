import meshpy.triangle as triangle
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse

"""
This module takes a boundary mesh and using the MeshPy package, creates
a mesh of the interior and immediate surroundings of that boundary. The
interior mesh is subdivided into regions depending on which side of a
boundary the cells are.
"""
class InteriorMesh(object):
    def __init__(self, facets, mesh_gen_scale_x_factor = 1.0, max_tri_area = None):
        meshpy_vs = facets.reshape((facets.shape[0] * facets.shape[1], facets.shape[2]))
        meshpy_facets = np.arange(meshpy_vs.shape[0]).reshape(
            (facets.shape[0], facets.shape[1])
        )
        meshpy_vs[:, 0] *= mesh_gen_scale_x_factor
        info = triangle.MeshInfo()
        info.set_points(meshpy_vs)
        info.set_facets(meshpy_facets)

        params = dict()
        if max_tri_area is not None:
            params['max_volume'] = max_tri_area
        mesh = triangle.build(info, **params)

        self.pts = np.array(mesh.points)
        self.pts[:, 0] /= mesh_gen_scale_x_factor
        self.tris = np.array(mesh.elements)
        self.boundary_facets = np.array(mesh.facets)

    def plot(self):
        plt.triplot(self.pts[:, 0], self.pts[:, 1], self.tris)
        plt.show()

    def region_plot(self):
        plt.tripcolor(
            self.pts[:, 0], self.pts[:, 1], self.tris,
            self.identify_regions()
        )
        plt.triplot(self.pts[:, 0], self.pts[:, 1], self.tris, 'r-')
        for f in self.boundary_facets:
            plt.plot(self.pts[f][:, 0], self.pts[f][:, 1], 'k')
        extents = expand_extents(
            *determine_extents(self.pts[self.boundary_facets]),
            factor = 1.1
        )
        plt.xlim(extents[0][0], extents[1][0])
        plt.ylim(extents[0][1], extents[1][1])
        plt.show()

    def identify_regions(self):
        n_tris = self.tris.shape[0]
        connectivity = scipy.sparse.dok_matrix((n_tris, n_tris))
        for t_idx in range(n_tris):
            adjacent_indices = find_adjacent_tris_indices(self.tris, self.tris[t_idx])
            for adj_idx in adjacent_indices:
                if self.are_across_an_input_facet(t_idx, adj_idx):
                    continue
                connectivity[t_idx, adj_idx] = 1

        n_components, components = scipy.sparse.csgraph.connected_components(
            connectivity, directed = False, return_labels = True
        )
        return components

    def are_across_an_input_facet(self, tri_A_idx, tri_B_idx):
        overlap = get_tri_pair_vertex_overlap(
            self.tris[tri_A_idx], self.tris[tri_B_idx]
        )
        pt_indices = [self.tris[tri_A_idx, o] for o in overlap]
        #TODO: Once this is too slow, build an index for fast searching
        for f in self.boundary_facets:
            if (f[0] == pt_indices[0] and f[1] == pt_indices[1]) or \
                (f[1] == pt_indices[0] and f[0] == pt_indices[1]):
                return True
        return False

def determine_extents(facets):
    min_corner = np.min(np.min(facets, axis = 0), axis = 0)
    max_corner = np.max(np.max(facets, axis = 0), axis = 0)
    return min_corner, max_corner

def expand_extents(min_corner, max_corner, factor = 2.0):
    center = (min_corner + max_corner) / 2.0
    dist = (max_corner - min_corner) / 2.0
    return center - dist * factor, center + dist * factor

def extents_to_box_2d(min_corner, max_corner):
    pts = [
        [min_corner[0], min_corner[1]],
        [max_corner[0], min_corner[1]],
        [max_corner[0], max_corner[1]],
        [min_corner[0], max_corner[1]]
    ]
    return np.array([
        [pts[0], pts[1]],
        [pts[1], pts[2]],
        [pts[2], pts[3]],
        [pts[3], pts[0]]
    ])

def add_extent_surface(input_facets):
    extents_box = extents_to_box_2d(
        *expand_extents(*determine_extents(input_facets))
    )
    return np.concatenate((input_facets, extents_box))


def find_adjacent_tris_indices(tris, query_tri):
    #TODO: Once this is too slow, build an index for fast searching
    return [i for i, t in enumerate(tris) if are_adjacent_tris(t, query_tri)]

def get_tri_pair_vertex_overlap(tri_A, tri_B):
    return [d for d in range(3) if tri_A[d] in tri_B]

def are_adjacent_tris(tri_A, tri_B):
    return len(get_tri_pair_vertex_overlap(tri_A, tri_B)) == 2
