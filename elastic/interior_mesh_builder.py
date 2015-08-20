import meshpy.triangle as triangle
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse

"""
This module takes a boundary mesh and using the MeshPy package, creates
a mesh of the interior and immediate surroundings of that boundary. The
interior mesh is subdivided into regions depending on which side of a
boundary the cells are.

The 'mesh_gen_scale_x_factor' allows the x dimension to be scaled in case
desired meshing should be well-proportioned in a different coordinate
scale than the input coordinates. For example, if the mesh will be used
for plotting and the x and y axes on the plot are not identical.
"""
class InteriorMeshBuilder(object):
    def __init__(self, facets, mesh_gen_scale_x_factor = 1.0, max_tri_area = None):
        meshpy_vs = facets.reshape((facets.shape[0] * facets.shape[1], facets.shape[2]))
        meshpy_facets = np.arange(meshpy_vs.shape[0]).reshape(
            (facets.shape[0], facets.shape[1])
        )

        meshpy_vs, meshpy_facets = self.remove_duplicate_vertices(
            meshpy_vs, meshpy_facets
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

    #TODO: This should probably already be done on the 3bem side, meshes
    # should be stored as vertex and cross-referenced triangle lists
    def remove_duplicate_vertices(self, meshpy_vs, meshpy_facets):
        equivalence_map = self.build_equivalence_map(meshpy_vs)
        new_vs = self.create_new_vertex_list(meshpy_vs, equivalence_map)
        new_facets = np.array(equivalence_map)[meshpy_facets]
        return new_vs, new_facets

    def build_equivalence_map(self, meshpy_vs):
        equivalence_map = []
        next = 0
        for v_idx1 in range(meshpy_vs.shape[0]):
            equivalence_map.append(next)
            next += 1
            for v_idx2 in range(v_idx1):
                sep_vec = meshpy_vs[v_idx1, :] - meshpy_vs[v_idx2, :]
                dist = np.sum(sep_vec ** 2)
                if np.allclose(meshpy_vs[v_idx1, :], meshpy_vs[v_idx2, :]):
                    equivalence_map[v_idx1] = equivalence_map[v_idx2]
                    next -= 1
                    break
        return equivalence_map

    def create_new_vertex_list(self, meshpy_vs, equivalence_map):
        out_vertices = np.max(equivalence_map) + 1
        new_vs = np.empty((out_vertices, 2))
        for v_idx1 in range(meshpy_vs.shape[0]):
            new_vs[equivalence_map[v_idx1], :] = meshpy_vs[v_idx1, :]
        return new_vs

    def plot(self, show = True):
        plt.triplot(self.pts[:, 0], self.pts[:, 1], self.tris)
        if show:
            plt.show()

    def region_plot(self, show = True):
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
        if show:
            plt.show()

    def identify_regions(self):
        #TODO: Once this is too slow, maybe build connectivity in c++ layer?
        n_tris = self.tris.shape[0]
        connectivity = scipy.sparse.dok_matrix((n_tris, n_tris))
        pts_to_tris_map = build_pt_to_tri_map(self.tris)
        for t_idx in range(n_tris):
            adjacent_indices = find_adjacent_tris_indices(
                self.tris, self.tris[t_idx], pts_to_tris_map
            )
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

def add_extent_surface(input_facets, expand_factor = 2.0):
    extents_box = extents_to_box_2d(
        *expand_extents(*determine_extents(input_facets), factor = expand_factor)
    )
    return np.concatenate((input_facets, extents_box))


def build_pt_to_tri_map(tris):
    pts_to_tris_map = dict()
    for t_idx, t in enumerate(tris):
        for d in range(3):
            pts_to_tris_map.setdefault(t[d], []).append(t_idx)
    return pts_to_tris_map

def find_adjacent_tris_indices(tris, query_tri, pts_to_tris_map):
    adjacent_tris = []
    for d in range(3):
        for touching_tri_idx in pts_to_tris_map[query_tri[d]]:
            if touching_tri_idx in adjacent_tris:
                continue
            if are_adjacent_tris(tris[touching_tri_idx], query_tri):
                adjacent_tris.append(touching_tri_idx)
    return adjacent_tris

def get_tri_pair_vertex_overlap(tri_A, tri_B):
    return [d for d in range(3) if tri_A[d] in tri_B]

def are_adjacent_tris(tri_A, tri_B):
    return len(get_tri_pair_vertex_overlap(tri_A, tri_B)) == 2
