import elastic

import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

def test_fault_inversion():
    fault_refine = 4
    es = []
    es.extend(elastic.line(
        [[-1, 0], [1, 0]], fault_refine,
        lambda pts: elastic.slip(pts, np.array([1 - np.abs(pts[:, 0]), [0, 0]]).T)
    ))

    params = dict(
        dense = True,
        src_far_order = 10
    )
    slip_imposed = elastic.execute(2, es, params)

    def get_normal_vec(pts):
        vector = np.array([
            -(pts[1][1] - pts[0][1]),
            pts[1][0] - pts[0][0]
        ])
        normal = vector / np.linalg.norm(vector)
        return normal

    def disp_bc(pts):
        disp = slip_imposed.interior_displacement(pts)
        return elastic.traction(pts, np.array(disp).T)

    def traction_bc(pts):
        normal = get_normal_vec(pts)
        trac = slip_imposed.interior_traction(pts, [normal, normal])
        return elastic.traction(pts, np.array(trac).T)

    def fault_bc(pts):
        normal = get_normal_vec(pts)
        trac = slip_imposed.interior_traction(pts, [normal, normal])
        return elastic.free_slip(pts, [0, 0], [trac[0][:]])

    bdry_refine = 5
    L = 4.0
    new_es = []
    new_es.extend(elastic.line([[-L, -L], [L, -L]], bdry_refine, disp_bc))
    new_es.extend(elastic.line([[L, -L], [L, L]], bdry_refine, traction_bc))
    new_es.extend(elastic.line([[L, L], [-L, L]], bdry_refine, traction_bc))
    new_es.extend(elastic.line([[-L, L], [-L, -L]], bdry_refine, traction_bc))
    crack_trac = np.array(slip_imposed.soln[('discontinuous', 'crack_traction')])
    fs = slip_imposed.meshes['discontinuous'].facets
    for i in range(fs.shape[0]):
        f = fs[i, :, :]
        pt_indices = [2 * i, 2 * i + 1]
        # new_es.append(fault_bc(f))
        new_es.append(elastic.crack(f, crack_trac[:, pt_indices].T))

    result = elastic.execute(2, new_es, params)

    slip_solved = np.array(result.soln[('discontinuous', 'slip')])
    slip_bc = np.array(slip_imposed.soln[('discontinuous', 'slip')])
    l2_error = np.sqrt(np.sum((slip_solved - slip_bc) ** 2))
    assert(l2_error / slip_solved.shape[1] < 0.002)
