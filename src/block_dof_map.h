#ifndef __1234jasHJAS_BLOCK_DOF_MAP_H
#define __1234jasHJAS_BLOCK_DOF_MAP_H

#include <map>

// Forward declarations
namespace tbem {
    template <size_t dim>
    struct Mesh;
}
template <size_t dim>
using MeshMap = std::map<std::string, tbem::Mesh<dim>>;

struct ComponentDescriptor
{
    const std::string where;
    const size_t component;
    bool operator<(const ComponentDescriptor& cd) const {
        if (where < cd.where) {
            return true;
        }
        if (where > cd.where) {
            return false;
        }
        return component < cd.component;
    }
};

struct BlockDOFMap {
    const std::map<ComponentDescriptor,size_t> map;
    const size_t total_dofs;
};

template <size_t dim>
BlockDOFMap from_mesh_set(const MeshMap<dim>& meshes);

#endif
