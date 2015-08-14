from elastic import displacement, traction, line
from elastic.interface import Executor
from elastic.mesh_provider import SkipUselessEntriesMeshProvider, SimpleMeshProvider

from elastic.bie_spec import get_elastic_kernels, get_BIEs
from elastic.compute import DenseIntegralEvaluator, IntegralDispatcher
from elastic.defaults import default_params
from elastic.dense_solver import DenseSolver
from elastic.iterative_solver import calculate_rhs

import tbempy

import numpy as np

def data(refine):
    es = []
    es.extend(
        line([[0, 0], [1, 0]], refine,
        lambda pts: traction(pts, [[1, 2], [3, 4]]))
    )
    es.extend(
        line([[1, 0], [2, 0]], refine,
        lambda pts: displacement(pts, [[1, 2], [3, 4]]))
    )
    executor = Executor(2, es, dict(dense = True))
    mp = SkipUselessEntriesMeshProvider(executor.meshes, executor.dof_map, [
        0, 1, 8, 9
    ])
    return mp, executor

def test_get_obs_mesh():
    mp, executor = data(1)
    mesh = mp.get_obs_mesh(dict(
        obs_mesh = 'continuous',
        unknown_field = 'traction'
    ))
    assert(mesh.n_facets() == 3)

def test_distribute_zeros():
    mp, executor = data(1)
    mesh = mp.get_obs_mesh(dict(
        obs_mesh = 'continuous',
        unknown_field = 'traction'
    ))
    A_np = np.random.rand(2 * mesh.n_dofs(), 2 * mesh.n_dofs())
    A = tbempy.TwoD.DenseOperator(A_np.shape[0], A_np.shape[1], A_np.reshape(A_np.size))
    B = mp.distribute_zeros(dict(
        obs_mesh = 'continuous',
        unknown_field = 'traction'
    ), A)
    B_np = B.data().reshape((B.n_rows(), B.n_cols()))
    for obs_dof in [0, 1, 8, 9]:
        np.testing.assert_almost_equal(B_np[obs_dof, :], 0)

def test_skip_with_matrix_construction():
    mp, executor = data(2)
    executor.mesh_provider = SimpleMeshProvider(executor.meshes)
    matrix_simple, rhs_simple = construct_final_system(executor)
    ignored_dofs = tbempy.TwoD.identify_ignored_dofs(executor.constraint_matrix)
    executor.mesh_provider = SkipUselessEntriesMeshProvider(
        executor.meshes, executor.dof_map, ignored_dofs
    )
    matrix_skipped, rhs_skipped = construct_final_system(executor)
    assert(np.all(matrix_simple - matrix_skipped == 0))
    assert(np.all(rhs_simple - rhs_skipped == 0))

def construct_final_system(executor):
    systems = executor.assemble(get_BIEs(executor.params))
    solver = DenseSolver(tbempy.TwoD, executor.params)
    matrix = solver.condense_matrix(
        solver.form_dense_matrix(executor.dof_map, systems),
        executor.constraint_matrix
    )
    rhs = calculate_rhs(
        tbempy.TwoD, executor.params, executor.dof_map,
        executor.constraint_matrix, systems
    )
    return matrix, rhs

