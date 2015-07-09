import numpy as np

class DOFMap(object):
    def __init__(self, dim, n_total_dofs, map):
        self.dim = dim
        self.n_total_dofs = n_total_dofs
        self.map = map

    @staticmethod
    def build(dim, field_types, meshes):
        map = dict()
        start = 0
        for mesh_name, field_name in field_types:
            components = []
            n_dofs = meshes[mesh_name].n_dofs()
            for d in range(dim):
                components.append(start)
                start += n_dofs
            components.append(start)
            map[(mesh_name, field_name)] = components
        n_total_dofs = start
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

        # displacement_field = fields[('continuous', 'displacement')]
        # fields[('continuous', 'gravity')] = [
        #     np.ones_like(displacement_field[d])
        #     for d in range(self.dim)
        # ]
        return fields

    def concatenate(self, fields):
        result = np.empty(self.n_total_dofs)
        for mesh_and_field, dofs in self.map.iteritems():
            for d in range(self.dim):
                result[dofs[d]:dofs[d+1]] = fields[mesh_and_field][d]
        return result
