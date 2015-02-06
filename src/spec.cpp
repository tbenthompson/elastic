#include "spec.h"
#include "filenames.h"
#include "3bem/continuity_builder.h"
#include "3bem/vertex_iterator.h"

using namespace tbem;

std::string trac_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".trac_out";
}

template <size_t dim>
ConstraintMatrix form_traction_constraints(const MeshMap<dim>& meshes, size_t d)
{
    return from_constraints({});
}

std::vector<IntegralSpec> displacement_BIE_terms(std::string obs_mesh) {
    IntegralSpec uut{obs_mesh, "displacement", "traction", "displacement", -1};
    IntegralSpec uuu{obs_mesh, "displacement", "displacement", "traction", -1};
    IntegralSpec utt{obs_mesh, "traction", "traction", "displacement", -1};
    IntegralSpec utu{obs_mesh, "traction", "displacement", "traction", -1};
    IntegralSpec ust{obs_mesh, "slip", "traction", "slip", -1};
    IntegralSpec uct{obs_mesh, "crack_traction", "traction", "slip", -1};
    // IntegralSpec uft{obs_mesh, "free_slip", "traction", "slip", -1};

    return {uut, utt, ust, uuu, utu, uct};
}

template <size_t dim>
IntegralEquationSpec<dim> get_displacement_BIE(const std::string& obs_mesh) {

    return {
        {obs_mesh, "displacement", 1},
        displacement_BIE_terms(obs_mesh),
        form_traction_constraints<dim>,
        trac_out_filename
    };
}

std::string disp_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".disp_out";
}

template <size_t dim>
ConstraintMatrix form_displacement_constraints(const MeshMap<dim>& meshes, size_t d)
{
    auto continuity = mesh_continuity(meshes.at("traction").begin());
    auto cut_continuity = cut_at_intersection(cut_at_intersection(
        continuity, meshes.at("traction").begin(), meshes.at("slip").begin()
    ), meshes.at("traction").begin(), meshes.at("crack_traction").begin());
    auto constraints = convert_to_constraints(cut_continuity);
    auto constraint_matrix = from_constraints(constraints);
    return constraint_matrix;
}

std::vector<IntegralSpec> traction_BIE_terms(std::string obs_mesh) {
    IntegralSpec tuh{obs_mesh, "displacement", "hypersingular", "displacement", 1};
    IntegralSpec tua{obs_mesh, "displacement", "adjoint_traction", "traction", 1};
    IntegralSpec tth{obs_mesh, "traction", "hypersingular", "displacement", 1};
    IntegralSpec tta{obs_mesh, "traction", "adjoint_traction", "traction", 1};
    IntegralSpec tsh{obs_mesh, "slip", "hypersingular", "slip", 1};
    IntegralSpec tch{obs_mesh, "crack_traction", "hypersingular", "slip", 1};
    // IntegralSpec tfh{obs_mesh, "free_slip", "hypersingular", "slip", 1};
    return {tuh, tth, tsh, tua, tta, tch};
}

template <size_t dim>
IntegralEquationSpec<dim> get_traction_BIE(const std::string& obs_mesh) {

    return {
        {obs_mesh, "traction", 1},
        traction_BIE_terms(obs_mesh),
        form_displacement_constraints<dim>,
        disp_out_filename
    };
}

template <size_t dim>
ConstraintMatrix form_slip_constraints(const MeshMap<dim>& meshes, size_t d)
{
    auto continuity = mesh_continuity(meshes.at("crack_traction").begin());
    auto constraints = convert_to_constraints(continuity);
    auto constraint_matrix = from_constraints(constraints);
    return constraint_matrix;
}

std::string slip_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".slip_out";
}

template <size_t dim>
IntegralEquationSpec<dim> get_crack_traction_BIE(const std::string& obs_mesh) {
    return {
        {obs_mesh, "crack_traction", 1},
        traction_BIE_terms(obs_mesh),
        form_slip_constraints<dim>,
        slip_out_filename
    };
}

std::vector<std::string> get_mesh_types() {
    return {
        "traction", "displacement", "slip", "crack_traction", "free_slip"
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
template 
IntegralEquationSpec<2> get_crack_traction_BIE(const std::string& obs_mesh);
template 
IntegralEquationSpec<3> get_crack_traction_BIE(const std::string& obs_mesh);
