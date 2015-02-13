from input_builder import Element, bem_template, run
import numpy as np

def main():
    vertices = np.array([
        [-1.0, -1.0,  1.0], [ 1.0, -1.0,  1.0], [ 1.0,  1.0,  1.0], [-1.0,  1.0,  1.0],
        [-1.0, -1.0, -1.0], [ 1.0, -1.0, -1.0], [ 1.0,  1.0, -1.0], [-1.0,  1.0, -1.0]
    ])
    facets = [
        # front
        [[0, 1, 2], [2, 3, 0], "traction", [0.0, 0.0, 0.0]],
        # top
        [[3, 2, 6], [6, 7, 3], "traction", [0.0, 0.0, 0.0]],
        # back
        [[7, 6, 5], [5, 4, 7], "traction", [0.0, 0.0, 0.0]],
        # bottom
        [[4, 5, 1], [1, 0, 4], "traction", [0.0, 0.0, 0.0]],
        # left
        [[4, 0, 3], [3, 7, 4], "displacement", [0.0, 0.0, 0.0]],
        # right
        [[1, 5, 6], [6, 2, 1], "displacement", [1.0, 0.0, 0.0]],
    ]

    es = []
    for f in facets:
        bc_type = f[2]
        bc = f[3]
        for t in [0, 1]:
            pts = [vertices[p_idx] for p_idx in f[t]]
            to_zero = -(pts[0] + pts[1] + pts[2]) / 3.0
            normal = np.cross(pts[0] - pts[1], pts[2] - pts[1])
            print np.sign(np.dot(normal, to_zero))
            es.append(Element([pts[1], pts[0], pts[2]], bc_type, [bc, bc, bc], 0))


    filename = "data/box3d.in"
    bem_template(filename, es = es, G = 1.0, mu = 0.25, solver_tol = 1e-6)
    run(filename, dim = 3)

if __name__ == '__main__':
    main()
