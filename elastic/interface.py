import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import input_builder
import compute
import dof_handling
import iterative_solver
import dense_solver

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

class Result(object):
    def __init__(self, tbem, soln, input):
        self.tbem = tbem
        self.soln = soln
        self.input = input

    """
    Compute the interior displacement at the specified points.
    """
    @log_exceptions
    def interior_displacement(self, pts):
        gravity = self.input.params['gravity']
        terms = bie_spec.displacement_BIE_terms("displacement", gravity)
        normals = np.zeros_like(pts)
        return compute.interior_eval(
            self.tbem, self.input.bcs, self.input.meshes, self.input.kernels,
            self.input.params, self.soln, pts, normals, terms
        )

    """
    Compute the interior traction on the planes specified by the normals and
    the specified points/
    """
    @log_exceptions
    def interior_traction(self, pts, normals):
        gravity = self.input.params['gravity']
        terms = bie_spec.traction_BIE_terms("traction", gravity)
        return compute.interior_eval(
            self.tbem, self.input.bcs, self.input.meshes, self.input.kernels,
            self.input.params, self.soln, pts, normals, terms
        )

    @log_exceptions
    def save(self, filename):
        with open(filename, 'w') as f:
            np.savez(
                f,
                dim = self.tbem.dim,
                soln = self.soln,
                elements = self.input.elements,
                params = self.input.params
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

@log_exceptions
def execute(dim, elements, input_params):
    data = assemble(dim, elements, input_params)
    return solve(*data)

def assemble(dim, elements, input_params):
    tbem = get_tbem(dim)
    input = input_builder.build_input(tbem, elements, input_params)
    dof_map = dof_handling.build_dof_map(tbem, input.bies, input.meshes)
    cs = dof_handling.build_constraint_matrix(tbem, dof_map, input)
    systems = compute.form_linear_systems(
        tbem, input.bies, input.meshes, input.bcs, input.kernels, input.params
    )
    return tbem, input, dof_map, cs, systems

def solve(tbem, input, dof_map, constraint_matrix, systems):
    solve_fnc = iterative_solver.iterative_solver
    if input.params['dense']:
        solve_fnc = dense_solver.dense_solver
    soln = solve_fnc(
        tbem, input.params, input.bies, dof_map, constraint_matrix, systems
    )
    return Result(tbem, soln, input)

def get_tbem(dim):
    if dim == 2:
        return tbempy.TwoD
    else:
        return tbempy.ThreeD
