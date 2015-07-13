import numpy as np

class DOFMap(object):
    def __init__(self, dim, n_total_dofs, map):
        self.dim = dim
        self.n_total_dofs = n_total_dofs
        self.map = map

    @staticmethod
    def build(dim, field_types, meshes):
        map = dict()
        next = 0
        for mesh_name, field_name in field_types:
            components = []
            n_dofs = meshes[mesh_name].n_dofs()
            for d in range(dim):
                components.append(next)
                next += n_dofs
            components.append(next)
            map[(mesh_name, field_name)] = components
        n_total_dofs = next
        return DOFMap(dim, n_total_dofs, map)

    def get(self, mesh_name, field_name, component_idx):
        return self.map[(mesh_name, field_name)][component_idx]

    def expand(self, distributed):
        fields = dict()
        for mesh_and_field, dofs in self.map.iteritems():
            field_components = []
            for d in range(self.dim):
                start_idx = dofs[d]
                end_idx = dofs[d + 1]
                field_components.append(distributed[start_idx:end_idx])
            fields[mesh_and_field] = field_components
        return fields

    def concatenate(self, fields):
        result = np.empty(self.n_total_dofs)
        for mesh_and_field, dofs in self.map.iteritems():
            for d in range(self.dim):
                if mesh_and_field in fields:
                    result[dofs[d]:dofs[d+1]] = fields[mesh_and_field][d]
                else:
                    result[dofs[d]:dofs[d+1]] = 0
        return result

    def get_matrix_block(self, output_type, input_type):
        row_components = self.map[output_type]
        col_components = self.map[input_type]
        start_row = row_components[0]
        end_row = row_components[-1]
        start_col = col_components[0]
        end_col = col_components[-1]
        return start_row, end_row, start_col, end_col
