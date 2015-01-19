#ifndef __reqrqweklajfbna_BEM_PARTS_H
#define __reqrqweklajfbna_BEM_PARTS_H
#include "load.h"
#include "spec.h"
#include "function.h"
#include "3bem/3bem.h"

template <size_t dim>
struct BEM {
    const Parameters params;
    const MeshMap<dim> meshes;
    const BCMap bcs;
    const KernelMap<dim> kernels;
    const tbem::QuadStrategy<dim> quad_strategy;
    const std::vector<IntegralEquationSpec> eqtn_specs;
    const std::vector<tbem::ConstraintMatrix> constraints;
};

struct ComputedOperator {
    const tbem::MatrixOperator op;
    const std::string src_mesh;
    const std::string function;
};

// The integral equation is implicitly defined by (sum(integrals) + mass) = 0
template <size_t dim>
std::vector<ComputedOperator>
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec& eqtn_spec);

struct LinearSystem {
    //TODO: Replace with ComputedIntegralEquation
    std::vector<ComputedOperator> lhs;
    Function rhs;
};

//TODO: better name
LinearSystem separate(const std::vector<ComputedOperator>& eqtn, const BCMap& bcs);

LinearSystem scale_rows(const LinearSystem& eqtn);

#endif
