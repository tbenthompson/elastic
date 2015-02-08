#include "UnitTest++.h"
#include "block_dof_map.h"
#include "3bem/mesh.h"

using namespace tbem;

TEST(FromMeshSet) {
    Facet<2> f{{{0.0, 0.0}, {0.0, 1.0}}};
    MeshMap<2> meshes; 
    meshes.insert(std::make_pair("def", Mesh<2>{{f,f,f}}));
    meshes.insert(std::make_pair("abc", Mesh<2>{{f}}));
    auto dof_map = from_mesh_set(meshes);
    CHECK_EQUAL((dof_map.map.at(ComponentDescriptor{"abc", 1})), 2);
    CHECK_EQUAL((dof_map.map.at(ComponentDescriptor{"def", 1})), 6 + 6);
}
