#ifndef __YYYCYYCXYZYZYYZI_SPEC_H
#define __YYYCYYCXYZYZYYZI_SPEC_H
#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <map>
#include "3bem/constraint_matrix.h"
#include "3bem/mesh.h"

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
using MeshMap = std::map<std::string, tbem::Mesh<dim>>;

template <size_t dim>
using ConstraintBuilder =
    std::function<tbem::ConstraintMatrix(const MeshMap<dim>&, size_t d)>;

// The convention will be that sum of the mass and the terms is equal to 0.
template <size_t dim>
struct IntegralEquationSpec {
    MassSpec mass;
    std::vector<IntegralSpec> terms;
    ConstraintBuilder<dim> constraint_builder;
};

template <size_t dim>
IntegralEquationSpec<dim> get_displacement_BIE(const std::string& obs_mesh);

template <size_t dim>
IntegralEquationSpec<dim> get_traction_BIE(const std::string& obs_mesh);

std::vector<std::string> get_mesh_types();
std::vector<std::string> get_bc_types();

#endif
