from elastic import line, traction, free_slip
from elastic.meshing import build_meshes
import tbempy

def test_split():
    for refine in range(7):
        es = []
        es.extend(line([[1, 1], [0, 1]], refine,
            lambda pts: traction(pts, [[0, 0], [0, 0]])
        ))
        n_continuous_facets = len(es)
        es.extend(line([[0.2, 0.3], [0.553182, 1.0]], refine,
            lambda pts: free_slip(pts, [0, 0], [[0, 0]])
        ))
        n_discontinuous_facets = len(es) - n_continuous_facets
        meshes = build_meshes(tbempy.TwoD, es)
        assert(meshes['continuous'].n_facets() == (n_continuous_facets + 1))
        assert(meshes['discontinuous'].n_facets() == n_discontinuous_facets)


if __name__ == '__main__':
    test_split()
