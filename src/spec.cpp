#include "spec.h"
#include "3bem/continuity_builder.h"
#include "3bem/vertex_iterator.h"

using namespace tbem;

template <size_t dim>
ConstraintMatrix form_traction_constraints(const MeshMap<dim>& meshes, size_t d)
{
    return from_constraints({});
}

template <size_t dim>
IntegralEquationSpec<dim> get_displacement_BIE(const std::string& obs_mesh) {
    IntegralSpec uut{obs_mesh, "displacement", "traction", "displacement", -1};
    IntegralSpec uuu{obs_mesh, "displacement", "displacement", "traction", -1};
    IntegralSpec utt{obs_mesh, "traction", "traction", "displacement", -1};
    IntegralSpec utu{obs_mesh, "traction", "displacement", "traction", -1};
    IntegralSpec ust{obs_mesh, "slip", "traction", "slip", -1};

    return {
        {obs_mesh, "displacement", 1},
        {uut, utt, ust, uuu, utu},
        form_traction_constraints<dim>
    };
}

template <size_t dim>
ConstraintMatrix form_displacement_constraints(const MeshMap<dim>& meshes, size_t d)
{
    auto continuity = mesh_continuity(meshes.at("traction").begin());
    auto cut_continuity = cut_at_intersection(
        continuity, meshes.at("traction").begin(), meshes.at("slip").begin()
    );
    auto constraints = convert_to_constraints(cut_continuity);
    auto constraint_matrix = from_constraints(constraints);
    return constraint_matrix;
}

template <size_t dim>
IntegralEquationSpec<dim> get_traction_BIE(const std::string& obs_mesh) {
    IntegralSpec tuh{obs_mesh, "displacement", "hypersingular", "displacement", 1};
    IntegralSpec tua{obs_mesh, "displacement", "adjoint_traction", "traction", 1};
    IntegralSpec tth{obs_mesh, "traction", "hypersingular", "displacement", 1};
    IntegralSpec tta{obs_mesh, "traction", "adjoint_traction", "traction", 1};
    IntegralSpec tsh{obs_mesh, "slip", "hypersingular", "slip", 1};

    return {
        {obs_mesh, "traction", 1},
        {tuh, tth, tsh, tua, tta},
        form_displacement_constraints<dim>
    };
}

std::vector<std::string> get_mesh_types() {
    return {
        "traction", "displacement", "slip"
    };
}

std::vector<std::string> get_bc_types() {
    return get_mesh_types();
}

template 
IntegralEquationSpec<2> get_displacement_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_displacement_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<2> get_traction_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_traction_BIE(const std::string& obs_mesh);
