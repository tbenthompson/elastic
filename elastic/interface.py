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
from exceptions import log_exceptions

import numpy as np
import logging
logger = logging.getLogger(__name__)


@log_exceptions
def execute(dim, elements, input_params):
    log_begin(len(elements), input_params)
    inputs = create_input(dim, elements, input_params)
    executor = Executor(*inputs)
    return executor.run()

def log_begin(n_elements, input_params):
    logger.info('')
    logger.info(
        'Beginning elastic calculation for ' +
        str(n_elements) +
        ' with params: ' +
        str(input_params)
    )

def create_input(dim, elements, input_params):
    arguments = (dim, elements, input_params)
    tbem = get_tbem(dim)
    params = defaults.add_default_parameters(input_params)
    meshes = meshing.build_meshes(tbem, elements)
    meshes = meshing.postprocess_meshes(tbem, meshes)
    return arguments, params, meshes

class Executor(object):
    def __init__(self, arguments, params, meshes):
        self.arguments = arguments
        self.params = params
        self.meshes = meshes
        self.setup()

    def setup(self):
        self.tbem = get_tbem(self.arguments[0])
        self.dof_map = dof_handling.DOFMap.build(
            self.tbem.dim, bie_spec.field_types, self.meshes
        )
        self.constraint_matrix = constraints.build_constraint_matrix(
            self.tbem, self.dof_map, self.arguments[1], self.meshes
        )

    def run(self):
        return self.solve(self.assemble(bie_spec.get_BIEs(self.params)))

    def assemble(self, bies):
        n_facets = str(self.meshes['all_mesh'].n_facets())
        logger.info('Assembling linear system for ' + n_facets + ' facets')

        evaluator = compute.FMMIntegralEvaluator(
            self.tbem, self.params, self.meshes['all_mesh']
        )
        if self.params['dense']:
            evaluator = compute.DenseIntegralEvaluator(
                self.tbem, self.params, self.meshes['all_mesh']
            )
        kernels = bie_spec.get_elastic_kernels(self.tbem, self.params)
        dispatcher = compute.IntegralDispatcher(self.meshes, kernels, evaluator)
        assembled_systems = system.form_linear_systems(bies, dispatcher)
        logger.info('Finished linear system assembly')
        return assembled_systems

    def solve(self, assembled_systems):
        solve_fnc = iterative_solver.iterative_solver
        if self.params['dense']:
            solve_fnc = dense_solver.dense_solver
        soln = solve_fnc(
            self.tbem, self.params, self.dof_map,
            self.constraint_matrix, assembled_systems
        )
        return Result(self.tbem, soln, self.meshes, self.params, self.arguments)

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
    @log_exceptions
    def interior_displacement(self, pts):
        normals = np.zeros_like(pts)
        return self._interior_eval(pts, normals, 'displacement')

    """
    Compute the interior traction at the points on the planes specified by the
    normals.
    """
    @log_exceptions
    def interior_traction(self, pts, normals):
        return self._interior_eval(pts, normals, 'traction')

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
        return system.evaluate_interior(dispatcher, self.soln, pts, normals, bie)

    @log_exceptions
    def save(self, filename):
        #TODO: at some point, include some information about versioning
        with open(filename, 'w') as f:
            np.savez(
                f,
                soln = self.soln,
                arguments = self.arguments
            )

    @staticmethod
    @log_exceptions
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
            meshes = meshing.build_meshes(tbem, elements)
        return Result(tbem, soln, meshes, params, arguments)

def get_tbem(dim):
    if dim == 2:
        return tbempy.TwoD
    else:
        return tbempy.ThreeD
