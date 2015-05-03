from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD
import bie_spec
import input_builder
from collections import namedtuple
import numpy as np

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def solve(dim, elements, input_params):
    input = input_builder.build_input(elements, input_params)
    tbem = get_tbem(dim)
    meshes, bcs = meshes_bcs_from_elements(tbem, elements)
    bies = bie_spec.get_all_BIEs()
    dof_map = build_block_dof_map([2, 3])

if __name__ == "__main__":
    solve(2, [
        Element(
            pts = [[0, 0], [1, 0]],
            bc = [[0, 0], [0, 0]],
            bc_type = 'displacement',
            n_refines = 0
        ),
        Element(
            pts = [[0, 0], [1, 0]],
            bc = [[1, 0], [1, 0]],
            bc_type = 'displacement',
            n_refines = 0
        )
    ], {})

