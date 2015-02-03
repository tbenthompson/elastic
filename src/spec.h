#ifndef __YYYCYYCXYZYZYYZI_SPEC_H
#define __YYYCYYCXYZYZYYZI_SPEC_H
#include <vector>
#include <string>
#include <iostream>

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

// The convention will be that sum of the mass and the terms is equal to 0.
struct IntegralEquationSpec {
    MassSpec mass;
    std::vector<IntegralSpec> terms;
};

IntegralEquationSpec get_displacement_BIE(const std::string& obs_mesh);

IntegralEquationSpec get_traction_BIE(const std::string& obs_mesh);

std::vector<std::string> get_mesh_types();
std::vector<std::string> get_bc_types();

#endif
