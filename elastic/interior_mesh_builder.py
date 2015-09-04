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
    def __init__(self, pts, tris, boundary_facets):
        self.pts = pts
        self.tris = tris
        self.boundary_facets = boundary_facets
        self.tri_region_map = tbempy.TwoD.identify_regions(
            self.tris.astype(np.uint64), self.boundary_facets.astype(np.uint64)
        )

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

    def refine(self, max_tri_area):
        #TODO: build_interior_mesh should take facets not mesh
        return build_interior_mesh(
            tbempy.TwoD.Mesh(self.pts[self.boundary_facets]),
            max_tri_area = max_tri_area
        )

    def get_region_id_from_pt(self, pt):
        result = tbempy.TwoD.find_containing_tri_idx(pt, self.tris, self.pts)
        if result.first is False:
            return None
        return self.tri_region_map[result.second]

    @log_elapsed_time(logger, 'retrieval of a specific subregion')
    def get_region(self, region_id):
        #TODO: Grab region via query point!
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
    tris = np.array(mesh.elements).astype(np.uint64)
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
