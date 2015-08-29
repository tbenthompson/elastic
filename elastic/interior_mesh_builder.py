from elastic.log_tools import log_elapsed_time
import tbempy.TwoD

import meshpy.triangle as triangle
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import time
import logging
logger = logging.getLogger(__name__)

"""
This module takes a boundary mesh and using the MeshPy package, creates
a mesh of the interior and immediate surroundings of that boundary. The
interior mesh is subdivided into regions depending on which side of a
boundary the cells are.

"""
class InteriorMesh(object):
    def __init__(self, pts, tris, boundary_facets, tri_region_map = None):
        self.pts = pts
        self.tris = tris
        self.boundary_facets = boundary_facets
        self.times = 0
        self.tri_region_map = tri_region_map
        if self.tri_region_map is None:
            self.tri_region_map = self.identify_regions()

    def plot(self, show = True):
        plt.triplot(self.pts[:, 0], self.pts[:, 1], self.tris)
        if show:
            plt.show()

    def region_plot(self, show = True):
        plt.tripcolor(
            self.pts[:, 0], self.pts[:, 1], self.tris,
            self.tri_region_map
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

    @log_elapsed_time(logger, 'identification of volumetric regions')
    def identify_regions(self):
        self.times += 1
        if self.times > 1:
            raise Exception()
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

    @log_elapsed_time(logger, 'retrieval of a specific subregion')
    def get_region(self, region_id):
        #TODO: Clean this up!
        #TODO: Grab multiple regions
        region_tris = self.tris[self.tri_region_map == region_id]
        out_pt_indices = np.unique(region_tris)

        tri_old_idx_to_new_idx_map = dict()
        for new_idx, old_idx in enumerate(out_pt_indices):
            tri_old_idx_to_new_idx_map[old_idx] = new_idx

        out_pts = np.empty((out_pt_indices.shape[0], 2))
        for old_idx in out_pt_indices:
            out_pts[tri_old_idx_to_new_idx_map[old_idx], :] = self.pts[old_idx]

        out_tris = np.empty(region_tris.shape)
        for i in range(region_tris.shape[0]):
            for d in range(region_tris.shape[1]):
                out_tris[i, d] = tri_old_idx_to_new_idx_map[region_tris[i, d]]

        out_boundary_facets = []
        for f_idx in range(self.boundary_facets.shape[0]):
            on_region_bdry = self.boundary_facets[f_idx, 0] in out_pt_indices and\
                self.boundary_facets[f_idx, 1] in out_pt_indices
            if on_region_bdry:
                out_facet = [
                    tri_old_idx_to_new_idx_map[old_idx] for old_idx in
                    self.boundary_facets[f_idx, :]
                ]
                out_boundary_facets.append(out_facet)

        return InteriorMesh(out_pts, out_tris, np.array(out_boundary_facets))


"""
The 'mesh_gen_scale_x_factor' allows the x dimension to be scaled in case
desired meshing should be well-proportioned in a different coordinate
scale than the input coordinates. For example, if the mesh will be used
for plotting and the x and y axes on the plot are not identical.
"""
@log_elapsed_time(logger, 'interior mesh construction')
def build_interior_mesh(bdry_mesh, mesh_gen_scale_x_factor = 1.0, max_tri_area = None):
    pt_index_mesh = tbempy.TwoD.convert_facet_to_pt_index(bdry_mesh)

    meshpy_vs = pt_index_mesh.points
    meshpy_facets = pt_index_mesh.facets
    meshpy_vs[:, 0] *= mesh_gen_scale_x_factor

    info = triangle.MeshInfo()
    info.set_points(meshpy_vs)
    info.set_facets(meshpy_facets.astype(np.int))

    params = dict()
    if max_tri_area is not None:
        params['max_volume'] = max_tri_area
    mesh = exec_triangle(info, params)

    pts = np.array(mesh.points)
    pts[:, 0] /= mesh_gen_scale_x_factor
    tris = np.array(mesh.elements)
    boundary_facets = np.array(mesh.facets)

    return InteriorMesh(pts, tris, boundary_facets)

@log_elapsed_time(logger, 'calling meshpy to construct interior mesh')
def exec_triangle(info, params):
    return triangle.build(info, **params)

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
