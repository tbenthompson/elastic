#ifndef __reqrqweklajfbna_BEM_PARTS_H
#define __reqrqweklajfbna_BEM_PARTS_H
#include "load.h"
#include "spec.h"
#include "function.h"
#include "3bem/3bem.h"

struct ComputedOperator {
    const tbem::MatrixOperator op;
    const std::string src_mesh;
    const std::string function;
};

// The integral equation is implicitly defined by (sum(integrals) + mass) = 0
struct ComputedIntegralEquation {
    const std::vector<ComputedOperator> terms;
};

template <size_t dim>
struct BEM {
    const MeshMap<dim> meshes;
    const BCMap bcs;
    const KernelMap<dim> kernels;
    const tbem::QuadStrategy<dim> quad_strategy;
    const tbem::ConstraintMatrix displacement_constraints;

    BEM(const Parameters& params, const MeshMap<dim>& meshes, const BCMap& bcs);

    static
    tbem::ConstraintMatrix form_displacement_constraints(const MeshMap<dim>& meshes);

    static 
    tbem::QuadStrategy<dim> form_quad_strategy(const Parameters& params);

    ComputedOperator compute_mass(const MassSpec& op_spec);

    ComputedOperator compute_integral(const IntegralSpec& op_spec);

    ComputedIntegralEquation
    compute_integral_equation(const IntegralEquationSpec& eqtn_spec);
};

struct LinearSystem {
    //TODO: Replace with ComputedIntegralEquation
    std::vector<ComputedOperator> lhs;
    Function rhs;
};

//TODO: better name
LinearSystem separate(const ComputedIntegralEquation& eqtn, const BCMap& bcs);

LinearSystem scale_rows(const LinearSystem& eqtn);

template <size_t dim>
BEM<dim> parse_into_bem(const std::string& filename);

#endif
