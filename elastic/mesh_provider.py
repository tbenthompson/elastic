from tbempy.TwoD import DenseOperator
import numpy as np

class SimpleMeshProvider(object):
    def __init__(self, meshes):
        self.meshes = meshes
    def get_obs_mesh(self, op_spec):
        return self.meshes[op_spec['obs_mesh']]
    def get_src_mesh(self, op_spec):
        return self.meshes[op_spec['src_mesh']]
    def distribute_zeros(self, op_spec, matrix):
        return matrix

class SkipUselessEntriesMeshProvider(object):
    def __init__(self, meshes, dof_map, ignored_dofs):
        self.meshes = meshes
        self.dim = self.meshes.values()[0].facets.shape[1]
        self.dof_map = dof_map
        self.ignored_dofs = ignored_dofs
        # self.constraint_matrix = constraint_matrix

    def get_obs_mesh(self, op_spec):
        full_mesh = self.meshes[op_spec['obs_mesh']]
        unnecessary_facets = []
        for facet_idx in range(full_mesh.n_facets()):
            if self.is_facet_constrained(op_spec, facet_idx):
                unnecessary_facets.append(facet_idx)
        restricted_mesh = full_mesh.remove_facets(unnecessary_facets)
        return restricted_mesh

    #TODO: unit tests for this code
    #TODO: Create some kind of op_spec structure that contains
    # functions for the input and output types
    #TODO: Seems like dof map should be in the c++ layer?
    def is_facet_constrained(self, op_spec, facet_idx):
        for basis_idx in range(self.dim):
            for field_dim in range(self.dim):
                offset = self.dof_map.get(
                    op_spec['obs_mesh'], op_spec['unknown_field'], field_dim
                )
                dof = offset + facet_idx * self.dim + basis_idx
                if dof not in self.ignored_dofs:
                    return False
        return True

    def distribute_zeros(self, op_spec, matrix):
        full_mesh = self.meshes[op_spec['obs_mesh']]
        new_matrix = np.empty((self.dim * full_mesh.n_dofs(), matrix.shape[1]))
        restricted_row = 0
        for field_dim in range(self.dim):
            component_offset = field_dim * full_mesh.n_dofs()
            for facet_idx in range(full_mesh.n_facets()):
                facet_offset = facet_idx * self.dim
                for basis_idx in range(self.dim):
                    obs_dof = component_offset + facet_offset + basis_idx
                    if self.is_facet_constrained(op_spec, facet_idx):
                        new_matrix[obs_dof, :] = 0
                    else:
                        new_matrix[obs_dof, :] = matrix[restricted_row, :]
                        restricted_row += 1
        return new_matrix

    def get_src_mesh(self, op_spec):
        return self.meshes[op_spec['src_mesh']]
