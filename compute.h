#ifndef __reqrqweklajfbna_BEM_PARTS_H
#define __reqrqweklajfbna_BEM_PARTS_H
#include "load.h"
#include "spec.h"
#include "function.h"
#include "system.h"
#include "3bem/operator.h"

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
ComputedIntegralEquation
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec& eqtn_spec);

struct LinearSystem {
    //TODO: Replace with ComputedIntegralEquation
    std::vector<ComputedOperator> lhs;
    Function rhs;
};

//TODO: better name
LinearSystem separate(const ComputedIntegralEquation& eqtn, const BCMap& bcs);

LinearSystem scale_rows(const LinearSystem& eqtn);

#endif
