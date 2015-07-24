import numpy as np
import tbempy

def build_meshes(tbem, es):
    element_lists = {k:[] for k in ['continuous', 'discontinuous']}
    for e in es:
        element_lists[e['type']].append(e)
    meshes = dict()
    for mesh_name,elements in element_lists.iteritems():
        pts = np.array([e['pts'] for e in elements])
        if len(elements) == 0:
            pts = np.empty((0, tbem.dim, tbem.dim))
        meshes[mesh_name] = tbem.Mesh(pts)
    return meshes, element_lists

def postprocess_meshes(tbem, meshes):
    meshes['continuous'] = split(tbem, meshes['continuous'], meshes['discontinuous'])
    meshes['discontinuous'] = split(tbem, meshes['discontinuous'], meshes['continuous'])
    meshes['all_mesh'] = tbem.Mesh.create_union(meshes.values())
    return meshes

def split(tbem, mesh_to_split, splitter):
    mp = tbem.MeshPreprocessor()
    intersections = mp.find_intersections(
        mesh_to_split.facets, splitter.facets
    )
    split_surface = mp.split_facets_at_intersections(
        mesh_to_split.facets, intersections
    )
    return tbem.Mesh(split_surface)

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
