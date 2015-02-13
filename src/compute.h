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
    const std::string obs_mesh;
    const std::string src_mesh;
    const std::string function;
};

typedef std::vector<ComputedOperator> ComputedEquation;

// The integral equation is implicitly defined by (sum(integrals) + mass) = 0
template <size_t dim>
ComputedEquation
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec<dim>& eqtn_spec);

struct NoDiagonalException {};

struct LinearSystem {
    const ComputedEquation lhs;
    const BlockFunction evaluated_terms;

    const ComputedOperator& get_diag_block() const {
        for (size_t lhs_idx = 0; lhs_idx < lhs.size(); lhs_idx++) {
            const auto& op = lhs[lhs_idx];
            if (op.src_mesh == op.obs_mesh) {
                return op;
            }
        }
        throw NoDiagonalException();
    }
};

LinearSystem 
evaluate_computable_terms(const ComputedEquation& eqtn, const FunctionMap& fields);

#endif
