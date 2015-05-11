import numpy as np

def l2_error(dim, mesh, exact_fnc, est_vec):
    l2_diff = np.zeros(dim)
    l2_exact = np.zeros(dim)
    for f_idx in range(mesh.facets.shape[0]):
        for v_idx in range(dim):
            global_v_idx = f_idx * dim + v_idx
            v = mesh.facets[f_idx, v_idx, :]
            exact = exact_fnc(v)
            for d_est in range(dim):
                est = est_vec[d_est][global_v_idx]
                l2_diff[d_est] += (exact[d_est] - est) ** 2
                l2_exact += exact[d_est] ** 2
    error = np.sqrt(l2_diff / l2_exact)
    return error

def check_error(problem, mesh_name, field_name, exact_fnc, limit):
    if problem.input.meshes[mesh_name].facets.shape[0] > 0:
        error = l2_error(problem.dim,
            problem.input.meshes[mesh_name],
            exact_fnc,
            problem.soln[(mesh_name, field_name)]
        )
        print('(mesh = ' + mesh_name + ', field = ' + field_name + ') error: '\
             + str(error))
        assert(np.all(error < limit))

