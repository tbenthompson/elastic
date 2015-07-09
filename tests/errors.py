import numpy as np
import logging

logger = logging.getLogger('elastic')

def l2_error(dim, pts, exact_fnc, est_vec):
    exact_example = exact_fnc([0] * dim)
    l2_diff = np.zeros(len(exact_example))
    l2_exact = np.zeros(len(exact_example))
    for v_idx in range(pts.shape[0]):
        v = pts[v_idx, :]
        exact = exact_fnc(v)
        for d_est in range(len(exact)):
            est = est_vec[d_est][v_idx]
            l2_diff[d_est] += (exact[d_est] - est) ** 2
            l2_exact += exact[d_est] ** 2
    error = np.sqrt(l2_diff / l2_exact)
    return error

def check_interior_error(pts, est, exact_fnc, limit):
    error = l2_error(pts.shape[1], pts, exact_fnc, est)
    logger.debug('Interior evaluation error: %s', error)
    assert(np.all(error < limit))

def check_error(result, mesh_name, field_name, exact_fnc, limit):
    if result.meshes[mesh_name].facets.shape[0] > 0:
        f = result.meshes[mesh_name].facets
        vs = f.reshape((f.shape[0] * f.shape[1], f.shape[2]))
        soln = result.soln[(mesh_name, field_name)]
        error = l2_error(result.tbem.dim, vs, exact_fnc, soln)
        logger.debug(
            'Boundary error for (mesh = %s, field = %s): %s',
             mesh_name, field_name, error
        )
        assert(np.all(error < limit))
