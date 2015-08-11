import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import compute
import dof_handling
import defaults
import constraints
import iterative_solver
import meshing
import dense_solver
import system
import mesh_provider
from log_tools import log_exceptions, log_elapsed_time

import numpy as np
import time
import logging
logger = logging.getLogger(__name__)


@log_exceptions(logger)
def execute(dim, elements, input_params):
    log_begin(len(elements), input_params)
    executor = Executor(dim, elements, input_params)
    return executor.run()

def log_begin(n_elements, input_params):
    logger.info('')
    logger.info(
        'Beginning elastic calculation for ' +
        str(n_elements) +
        ' with params: ' +
        str(input_params)
    )

class Executor(object):
    def __init__(self, *data):
        self.set_params(*data)
        self.setup(*data)

    def set_params(self, dim, elements, input_params):
        self.arguments = (dim, elements, input_params)
        params = defaults.add_default_parameters(input_params)
        self.check_input_params(params)
        self.params = params

    @log_elapsed_time(logger, 'setup of BEM problem')
    def setup(self, dim, elements, input_params):
        self.tbem = get_tbem(dim)
        self.meshes, self.element_lists = meshing.build_meshes(self.tbem, elements)
        self.dof_map = dof_handling.DOFMap.build(
            self.tbem.dim, bie_spec.field_types, self.meshes
        )
        self.constraint_matrix = constraints.build_constraint_matrix(
            self.tbem, self.dof_map, self.arguments[1], self.meshes
        )
        ignored_dofs = self.tbem.identify_ignored_dofs(self.constraint_matrix)
        self.mesh_provider = mesh_provider.SkipUselessEntriesMeshProvider(
            self.meshes, self.dof_map, ignored_dofs
        )

    def check_input_params(self, params):
        if 'obs_order' in params:
            raise Exception(
                'obs_order parameter has been removed in ' +
                'favor of obs_near_order and obs_far_order'
            )
        if params['obs_far_order'] == params['src_far_order']:
            raise Exception(
                'obs_far_order and src_far_order can\'t be equal. ' +
                ' Please change one.'
            )

    def run(self):
        return self.solve(self.assemble(bie_spec.get_BIEs(self.params)))

    @log_elapsed_time(logger, 'assembly of BEM problem')
    def assemble(self, bies):
        evaluator = compute.FMMIntegralEvaluator(
            self.tbem, self.params, self.meshes['all_mesh']
        )
        if self.params['dense']:
            evaluator = compute.DenseIntegralEvaluator(
                self.tbem, self.params, self.meshes['all_mesh']
            )
        kernels = bie_spec.get_elastic_kernels(self.tbem, self.params)

        dispatcher = compute.IntegralDispatcher(self.mesh_provider, kernels, evaluator)
        assembled_systems = system.form_linear_systems(bies, dispatcher)
        return assembled_systems

    @log_elapsed_time(logger, 'linear system solution of BEM problem')
    def solve(self, assembled_systems):
        if self.params['dense']:
            solver = dense_solver.DenseSolver(self.tbem, self.params)
            soln = solver.solve(self.dof_map, self.constraint_matrix, assembled_systems)
        else:
            soln = iterative_solver.iterative_solver(
                self.tbem, self.params, self.dof_map,
                self.constraint_matrix, assembled_systems
            )
        out = Result(self.tbem, soln, self.meshes, self.params, self.arguments)
        return out

class Result(object):
    def __init__(self, tbem, soln, meshes, params, arguments):
        self.tbem = tbem
        self.soln = soln
        self.meshes = meshes
        self.params = params
        self.arguments = arguments

    """
    Compute the interior displacement at the specified points.
    """
    @log_exceptions(logger)
    def interior_displacement(self, pts):
        normals = np.zeros_like(pts)
        return self._interior_eval(pts, normals, 'displacement')

    """
    Compute the interior traction at the points on the planes specified by the
    normals.
    """
    @log_exceptions(logger)
    def interior_traction(self, pts, normals):
        return self._interior_eval(pts, normals, 'traction')

    @log_elapsed_time(logger, lambda self, args: 'interior evaluation of '
        + str(args[2]) + 's at ' + str(len(args[0])) + ' points')
    def _interior_eval(self, pts, normals, which_bie):
        gravity = self.params['gravity']
        bie = bie_spec.bie_from_field_name(
            'continuous', bie_spec.unknowns_to_knowns[which_bie], self.params
        )
        evaluator = compute.FMMIntegralEvaluator(
            self.tbem, self.params, self.meshes['all_mesh']
        )
        if self.params['dense']:
            evaluator = compute.DenseIntegralEvaluator(
                self.tbem, self.params, self.meshes['all_mesh']
            )
        kernels = bie_spec.get_elastic_kernels(self.tbem, self.params)
        dispatcher = compute.IntegralDispatcher(self.meshes, kernels, evaluator)
        system.add_constant_fields(self.soln, False)
        out = system.evaluate_interior(dispatcher, self.soln, pts, normals, bie)
        return out

    @log_exceptions(logger)
    def save(self, filename):
        #TODO: at some point, include some information about versioning
        with open(filename, 'w') as f:
            np.savez(
                f,
                soln = self.soln,
                arguments = self.arguments
            )

    @staticmethod
    @log_exceptions(logger)
    def load(filename):
        with open(filename, 'r') as f:
            npzfile = np.load(f)
            arguments = npzfile['arguments']
            dim = arguments[0]
            elements = arguments[1]
            input_params = arguments[2]
            soln = npzfile['soln'].tolist()
            tbem = get_tbem(dim)
            params = defaults.add_default_parameters(input_params)
            meshes, _ = meshing.build_meshes(tbem, elements)
        return Result(tbem, soln, meshes, params, arguments)

def get_tbem(dim):
    if dim == 2:
        return tbempy.TwoD
    else:
        return tbempy.ThreeD
