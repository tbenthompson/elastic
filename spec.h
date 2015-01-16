#ifndef __YYYCYYCXYZYZYYZI_SPEC_H
#define __YYYCYYCXYZYZYYZI_SPEC_H

#include <vector>
#include <string>
#include <memory>
#include <map>
#include "3bem/kernel.h"

template <size_t dim>
using KernelPtr = 
    std::shared_ptr<tbem::Kernel<dim,
        tbem::Vec<double,dim>,
        tbem::Vec<double,dim>,
        tbem::Vec<tbem::Vec<double,dim>,dim>
    >>;

template <size_t dim>
using KernelSet = std::map<std::string, KernelPtr<dim>>;

template <size_t dim>
KernelSet<dim> get_elastic_kernels(double shear_modulus, double poisson_ratio);

struct OperatorSpec 
{
    const std::string obs_mesh;
    const std::string src_mesh;
    const std::string kernel;
    const std::string function;
    const double multiplier;

    friend std::ostream& operator<<(std::ostream& os, const OperatorSpec& spec) {
        os << "(" << spec.obs_mesh << ", " << spec.src_mesh << ", " <<
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
    std::vector<OperatorSpec> terms;
};

IntegralEquationSpec get_displacement_BIE();

IntegralEquationSpec get_traction_BIE();

#endif
