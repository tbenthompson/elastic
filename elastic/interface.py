import numpy as np

import tbempy.TwoD
import tbempy.ThreeD

import bie_spec
import input_builder
import compute
import dof_handling
import iterative_solver
import dense_solver


class Result(object):
    def __init__(self, tbem, soln, input):
        self.tbem = tbem
        self.soln = soln
        self.input = input

    """
    Compute the interior displacement at the specified points.
    """
    def interior_displacement(self, pts):
        bie = bie_spec.get_displacement_BIE("displacement", self.input.params)
        normals = np.zeros_like(pts)
        return compute.interior_eval(
            self.tbem, self.input, self.soln, pts, normals, bie
        )

    """
    Compute the interior traction on the planes specified by the normals and
    the specified points/
    """
    def interior_traction(self, pts, normals):
        bie = bie_spec.get_traction_BIE("traction", self.input.params)
        return compute.interior_eval(
            self.tbem, self.input, self.soln, pts, normals, bie
        )

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

def execute(dim, elements, input_params):
    data = assemble(dim, elements, input_params)
    return solve(*data)

def assemble(dim, elements, input_params):
    dim = dim
    elements = elements
    input_params = input_params
    tbem = get_tbem(dim)
    input = input_builder.build_input(tbem, elements, input_params)
    dof_map = dof_handling.build_dof_map(tbem, input.bies, input.meshes)
    constraint_matrix = dof_handling.build_constraint_matrix(
        tbem, dof_map, input
    )
    systems = compute.form_linear_systems(tbem, input)
    return tbem, input, dof_map, constraint_matrix, systems

def solve(tbem, input, dof_map, constraint_matrix, systems):
    solve_fnc = iterative_solver.iterative_solver
    if input.params['dense']:
        solve_fnc = dense_solver.dense_solver
    soln = solve_fnc(
        tbem, input, dof_map, constraint_matrix, systems
    )
    return Result(tbem, soln, input)

def get_tbem(dim):
    if dim == 2:
        return tbempy.TwoD
    else:
        return tbempy.ThreeD
