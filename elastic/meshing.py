import numpy as np
import tbempy

def build_meshes(tbem, mesh_types, es):
    element_lists = {k:[] for k in mesh_types}
    for e in es:
        element_lists[e['type']].append(e)
    meshes = dict()
    for mesh_name,elements in element_lists.iteritems():
        pts = np.array([e['pts'] for e in elements])
        if len(elements) == 0:
            pts = np.empty((0, tbem.dim, tbem.dim))
        meshes[mesh_name] = tbem.Mesh(pts)

    mp = tbem.MeshPreprocessor()
    intersections = mp.find_intersections(
        meshes['continuous'].facets, meshes['discontinuous'].facets
    )
    split_surface = mp.split_facets_at_intersections(
        meshes['continuous'].facets, intersections
    )
    meshes['continuous'] = tbem.Mesh(split_surface)
    meshes['all_mesh'] = tbem.Mesh.create_union(meshes.values())
    return meshes

def line(end_pts, refine, element_builder):
    mesh = tbempy.TwoD.line_mesh(end_pts[0], end_pts[1]).refine_repeatedly(refine)
    return add_bcs(mesh, element_builder)

def circle(center, r, refine, element_builder, reverse):
    mesh = tbempy.TwoD.circle_mesh(center, r, refine, reverse)
    return add_bcs(mesh, element_builder)

def sphere(center, r, refine, element_builder, reverse):
    mesh = tbempy.ThreeD.sphere_mesh(center, r, refine, reverse)
    return add_bcs(mesh, element_builder)

def add_bcs(mesh, element_builder):
    return [
        element_builder(np.array(
            [mesh.facets[i, d, :] for d in range(mesh.facets.shape[1])]
        ))
        for i in range(mesh.n_facets())
    ]
