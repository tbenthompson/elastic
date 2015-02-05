#ifndef __reqrqweklajfbna_BEM_PARTS_H
#define __reqrqweklajfbna_BEM_PARTS_H
#include "load.h"
#include "kernels.h"
#include "3bem/operator.h"

template <size_t dim>
struct BEM;
template <size_t dim>
struct IntegralEquationSpec;

struct ComputedOperator {
    const tbem::BlockOperator op;
    const std::string src_mesh;
    const std::string function;
};

typedef std::vector<ComputedOperator> ComputedIntegralEquation;

// The integral equation is implicitly defined by (sum(integrals) + mass) = 0
template <size_t dim>
ComputedIntegralEquation
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec<dim>& eqtn_spec);

struct LinearSystem {
    ComputedIntegralEquation lhs;
    BlockFunction rhs;
};

//TODO: better name
LinearSystem separate(const ComputedIntegralEquation& eqtn, const BCMap& bcs);

LinearSystem scale_rows(const LinearSystem& eqtn);

#endif
