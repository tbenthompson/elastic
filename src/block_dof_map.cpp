#include "block_dof_map.h"
#include "load.h"

template <size_t dim>
BlockDOFMap from_mesh_set(const MeshMap<dim>& meshes) {
    std::map<ComponentDescriptor,size_t> dof_map;
    size_t n_cumulative_dofs = 0;
    for (auto it = meshes.begin(); it != meshes.end(); ++it) {
        const auto& m = it->second;
        auto n_dofs = m.n_dofs();
        for (size_t d = 0; d < 3; d++) {
            ComponentDescriptor component{it->first, d};
            dof_map[component] = n_cumulative_dofs;
            n_cumulative_dofs += n_dofs; 
        }
    }
    return {dof_map, n_cumulative_dofs};
}

template 
BlockDOFMap from_mesh_set(const MeshMap<2>& meshes);
template 
BlockDOFMap from_mesh_set(const MeshMap<3>& meshes);
