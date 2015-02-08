#ifndef __YYYCYYCXYZYZYYZI_SPEC_H
#define __YYYCYYCXYZYZYYZI_SPEC_H
#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <map>
#include "block_dof_map.h"
#include "3bem/constraint_matrix.h"
#include "3bem/mesh.h"

//Forward declarations
template <size_t dim>
using MeshMap = std::map<std::string, tbem::Mesh<dim>>;

struct IntegralSpec 
{
    const std::string obs_mesh;
    const std::string src_mesh;
    const std::string kernel;
    const std::string function;
    const double multiplier;

    friend std::ostream& operator<<(std::ostream& os, const IntegralSpec& spec) {
        os << "IntegralSpec(" << spec.obs_mesh << ", " << spec.src_mesh << ", " <<
            spec.kernel << ", " << spec.function << ", " <<  spec.multiplier << ")";
        return os;
    }
};

struct MassSpec {
    const std::string obs_mesh;
    const std::string function;
    const double multiplier;
};

template <size_t dim>
using ConstraintBuilder = std::function<std::vector<tbem::ConstraintEQ>
    (const MeshMap<dim>&, size_t d)>;

// The convention will be that sum of the mass and the terms is equal to 0.
template <size_t dim>
struct IntegralEquationSpec {
    std::string unknown_field;
    MassSpec mass;
    std::vector<IntegralSpec> terms;
    ConstraintBuilder<dim> constraint_builder;
    std::function<std::string(const std::string&)> get_output_filename;

    std::string obs_mesh() const {
        return mass.obs_mesh;
    }
};

template <size_t dim>
IntegralEquationSpec<dim> get_displacement_BIE(const std::string& obs_mesh);

template <size_t dim>
IntegralEquationSpec<dim> get_traction_BIE(const std::string& obs_mesh);

template <size_t dim>
IntegralEquationSpec<dim> get_crack_traction_BIE(const std::string& obs_mesh);

template <size_t dim>
std::vector<IntegralEquationSpec<dim>> get_all_BIEs() {
    return {
        get_displacement_BIE<dim>("displacement"),
        get_traction_BIE<dim>("traction"),
        get_crack_traction_BIE<dim>("crack_traction")
    };
}

std::vector<std::string> get_mesh_types();
std::vector<std::string> get_bc_types();

#endif
