import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import compute
import dof_handling
import defaults
import constraints
import iterative_solver
import meshing
# import dense_solver

import numpy as np
import traceback
import logging
logger = logging.getLogger(__name__)

def log_exceptions(f):
    def _wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            upper = ''.join(traceback.format_list(traceback.extract_stack()[:-1]))
            lower = '\n'.join(traceback.format_exc().split('\n')[1:-1])
            logger.error('Exception:\n' + upper + lower)
            raise
    return _wrapper

@log_exceptions
def execute(dim, elements, input_params):
    executor = Executor(dim, elements, input_params)
    return executor.solve()

class Executor(object):
    def __init__(self, dim, elements, input_params):
        self.arguments = (dim, elements, input_params)
        self.tbem = get_tbem(dim)
        self.params = defaults.add_default_parameters(input_params)
        self.meshes = meshing.build_meshes(
            self.tbem, bie_spec.mesh_types, elements
        )
        self.dof_map = dof_handling.DOFMap.build(
            self.tbem.dim, bie_spec.field_types, self.meshes
        )
        self.constraint_matrix = constraints.build_constraint_matrix(
            self.tbem, self.dof_map, elements, self.meshes
        )
        self.bies = bie_spec.get_BIEs(self.params)
        self.assemble()

    def assemble(self):
        evaluator = compute.DenseIntegralEvaluator(
            self.tbem, self.params, self.meshes['all_mesh']
        )
        kernels = bie_spec.get_elastic_kernels(self.tbem, self.params)
        dispatcher = compute.IntegralDispatcher(self.meshes, kernels, evaluator)
        self.systems = compute.form_linear_systems(self.bies, dispatcher)

    def solve(self):
        solve_fnc = iterative_solver.iterative_solver
        if self.params['dense']:
            solve_fnc = dense_solver.dense_solver
        soln = solve_fnc(
            self.tbem, self.params, self.bies, self.dof_map,
            self.constraint_matrix, self.systems
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
        gravity = self.params['gravity']
        terms = bie_spec.displacement_BIE_terms("displacement", gravity)
        normals = np.zeros_like(pts)
        kernels = bie_spec.get_elastic_kernels(self.tbem, self.params)
        return compute.interior_eval(
            self.tbem, self.meshes, kernels,
            self.params, self.soln, pts, normals, terms
        )

    """
    Compute the interior traction on the planes specified by the normals and
    the specified points/
    """
    @log_exceptions
    def interior_traction(self, pts, normals):
        return self.interior_eval(pts, normals, 'traction')

    def _interior_eval(self, pts, normals, which_bie):
        gravity = self.params['gravity']
        terms = bie_spec.traction_BIE_terms("traction", gravity)
        kernels = bie_spec.get_elastic_kernels(self.tbem, self.params)
        return compute.interior_eval(
            self.tbem, self.meshes, self.kernels,
            self.params, self.soln, pts, normals, terms
        )

    @log_exceptions
    def save(self, filename):
        with open(filename, 'w') as f:
            np.savez(
                f,
                dim = self.tbem.dim,
                soln = self.soln,
                elements = self.elements,
                params = self.params
            )

    @staticmethod
    @log_exceptions
    def load(filename):
        with open(filename, 'r') as f:
            npzfile = np.load(f)
            soln = npzfile['soln'].tolist()
            raw_elements = npzfile['elements'].tolist()
            elements = [
                input_builder.Element(*r) for r in raw_elements
            ]
            params = npzfile['params'].tolist()
            dim = npzfile['dim']
        tbem = get_tbem(dim)
        input = input_builder.build_input(tbem, elements, params)
        return Result(tbem, soln, input)

def get_tbem(dim):
    if dim == 2:
        return tbempy.TwoD
    else:
        return tbempy.ThreeD
