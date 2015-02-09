#include "spec.h"
#include "filenames.h"
#include "3bem/constraint_builder.h"
#include "3bem/vertex_iterator.h"
#include "3bem/block_dof_map.h"

using namespace tbem;

std::vector<ConstraintEQ> 
shift_constraints(const std::vector<ConstraintEQ>& constraints, size_t shift_dof) {
    std::vector<ConstraintEQ> out;
    for (size_t i = 0; i < constraints.size(); i++) {
        auto n_terms = constraints[i].terms.size();
        std::vector<LinearTerm> out_terms; 
        for (size_t t_idx = 0; t_idx < n_terms; t_idx++) {
            out_terms.push_back({
                constraints[i].terms[t_idx].dof + shift_dof,
                constraints[i].terms[t_idx].weight
            });
        }
        out.push_back({out_terms, constraints[i].rhs});
    }

    return out;
}

std::vector<std::string> get_mesh_types() {
    return {
        "traction", "displacement", "slip", "crack_traction", "free_slip_traction"
    };
}

std::vector<std::string> get_bc_types() {
    return get_mesh_types();
}

// ---- Displacement BIE

std::string trac_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".trac_out";
}


template <size_t dim>
std::vector<ConstraintEQ> form_traction_constraints(
    const std::map<std::string,size_t>& component_map,
    const BlockDOFMap& dof_map,
    const MeshMap<dim>& meshes,
    size_t d)
{
    return {};
}


std::vector<IntegralSpec> displacement_BIE_terms(std::string obs_mesh) {
    return {
        {obs_mesh, "displacement", "traction", "displacement", -1},
        {obs_mesh, "displacement", "displacement", "traction", -1},
        {obs_mesh, "traction", "traction", "displacement", -1},
        {obs_mesh, "traction", "displacement", "traction", -1},
        {obs_mesh, "slip", "traction", "slip", -1},
        {obs_mesh, "crack_traction", "traction", "slip", -1},
        {obs_mesh, "free_slip_traction", "traction", "free_slip", -1}
    };
}

template <size_t dim>
IntegralEquationSpec<dim> get_displacement_BIE(const std::string& obs_mesh) {

    return {
        "traction",
        {obs_mesh, "displacement", 1},
        displacement_BIE_terms(obs_mesh),
        form_traction_constraints<dim>,
        trac_out_filename
    };
}

// ---- Traction BIE

std::string disp_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".disp_out";
}

template <size_t dim>
std::vector<ConstraintEQ> form_displacement_constraints(
        const std::map<std::string,size_t>& component_map,
        const BlockDOFMap& dof_map,
        const MeshMap<dim>& meshes,
        size_t d)
{
    auto continuity = mesh_continuity(meshes.at("traction").begin());
    auto cut_continuity = cut_at_intersection(cut_at_intersection(
        continuity, meshes.at("traction").begin(), meshes.at("slip").begin()
    ), meshes.at("traction").begin(), meshes.at("crack_traction").begin());
    auto constraints = convert_to_constraints(cut_continuity);
    auto shift_dof = dof_map.start_positions[component_map.at("traction") + d];
    return shift_constraints(constraints, shift_dof);
}

std::vector<IntegralSpec> traction_BIE_terms(std::string obs_mesh) {
    return {
        {obs_mesh, "displacement", "hypersingular", "displacement", 1},
        {obs_mesh, "displacement", "adjoint_traction", "traction", 1},
        {obs_mesh, "traction", "hypersingular", "displacement", 1},
        {obs_mesh, "traction", "adjoint_traction", "traction", 1},
        {obs_mesh, "slip", "hypersingular", "slip", 1},
        {obs_mesh, "crack_traction", "hypersingular", "slip", 1},
        {obs_mesh, "free_slip_traction", "hypersingular", "free_slip", 1}
    };
}

template <size_t dim>
IntegralEquationSpec<dim> get_traction_BIE(const std::string& obs_mesh) {

    return {
        "displacement",
        {obs_mesh, "traction", 1},
        traction_BIE_terms(obs_mesh),
        form_displacement_constraints<dim>,
        disp_out_filename
    };
}

// ---- Crack traction BIE

template <size_t dim>
std::vector<ConstraintEQ> form_slip_constraints(
    const std::map<std::string,size_t>& component_map,
    const BlockDOFMap& dof_map,
    const MeshMap<dim>& meshes,
    size_t d)
{
    auto continuity = mesh_continuity(meshes.at("crack_traction").begin());
    auto constraints = convert_to_constraints(continuity);
    auto shift_dof = dof_map.start_positions[component_map.at("crack_traction") + d];
    return shift_constraints(constraints, shift_dof);
}

std::string slip_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".slip_out";
}

template <size_t dim>
IntegralEquationSpec<dim> get_crack_traction_BIE(const std::string& obs_mesh) {
    return {
        "slip",
        {obs_mesh, "crack_traction", 1},
        traction_BIE_terms(obs_mesh),
        form_slip_constraints<dim>,
        slip_out_filename
    };
}

// ---- Free slip crack traction BIE

template <size_t dim>
std::vector<ConstraintEQ> form_free_slip_constraints(
    const std::map<std::string,size_t>& component_map,
    const BlockDOFMap& dof_map,
    const MeshMap<dim>& meshes,
    size_t d)
{
    if (d != 0) {
        return {};
    }

    const auto& mesh = meshes.at("free_slip_traction");
    auto fs_component = component_map.at("free_slip_traction");
    auto x_start = dof_map.start_positions[fs_component];
    auto y_start = dof_map.start_positions[fs_component + 1];
    // Opening displacement constraint.
    std::vector<ConstraintEQ> constraints;
    for (auto it = mesh.begin(); it != mesh.end(); ++it) {
        auto ux_dof = x_start + it.absolute_index();
        auto uy_dof = y_start + it.absolute_index();
        auto normal = normalized(unscaled_normal(it.get_facet()));
        ConstraintEQ c{{{ux_dof, normal[0]}, {uy_dof, normal[1]}}, 0.0};
        constraints.push_back(c);
    }
    return constraints;
}

std::string free_slip_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".free_slip_out";
}

template <size_t dim>
IntegralEquationSpec<dim> get_free_slip_BIE(const std::string& obs_mesh) {
    return {
        "free_slip",
        {obs_mesh, "free_slip_traction", 1},
        traction_BIE_terms(obs_mesh),
        form_free_slip_constraints<dim>,
        free_slip_out_filename
    };
}

template 
IntegralEquationSpec<2> get_displacement_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_displacement_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<2> get_traction_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_traction_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<2> get_crack_traction_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_crack_traction_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<2> get_free_slip_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_free_slip_BIE(const std::string& obs_mesh);
