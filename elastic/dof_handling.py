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
        for mesh_name, field_name in field_types.keys():
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

        displacement_field = fields[('continuous', 'displacement')]
        return fields

    def concatenate(self, fncs):
        return np.concatenate(fncs)

#TODO: These don't really make sense in this location, where should they go? To an "evaluate" module?
def scale_columns(unknowns, scaling_fncs, params):
    for u, values in unknowns.iteritems():
        f = scaling_fncs[u]
        for d in range(len(values)):
            values[d] *= f(params)

def scale_rows(eval, bies, params):
    for i, bie in enumerate(bies):
        eval[i] *= bie['scaling'](params)
