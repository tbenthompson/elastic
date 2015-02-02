#ifndef __reqrqweklajfbna_BEM_PARTS_H
#define __reqrqweklajfbna_BEM_PARTS_H
#include "load.h"
#include "spec.h"
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
    const tbem::BlockOperator op;
    const std::string src_mesh;
    const std::string function;
};

// The integral equation is implicitly defined by (sum(integrals) + mass) = 0
template <size_t dim>
std::vector<ComputedOperator>
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec& eqtn_spec);

template <size_t dim>
std::vector<std::vector<double>>
compute_interior(const std::vector<tbem::ObsPt<dim>>& pts,
                 const BEM<dim>& bem,
                 const IntegralEquationSpec& eqtn_spec,
                 const BCMap& bcs) {
    std::vector<std::vector<double>> results(dim, std::vector<double>(pts.size(), 0.0));
    for (size_t i = 0; i < pts.size(); i++) {
        for (const auto& term: eqtn_spec.terms) {
            const auto& src_mesh = bem.meshes.at(term.src_mesh);
            const auto& kernel = bem.kernels.at(term.kernel);
            const auto& field =
                bcs.at(FieldDescriptor{term.src_mesh, term.function});

            auto problem = make_problem(src_mesh, tbem::Mesh<dim>{{}}, *kernel);
            auto op = mesh_to_point_operator(problem, bem.quad_strategy, pts[i]);
            auto evaluated = apply_operator(op, field);
            for (size_t d = 0; d < dim; d++) {
                results[d][i] += evaluated[d][0];
            }
        }
    }
    return results;
}

struct LinearSystem {
    //TODO: Replace with ComputedIntegralEquation
    std::vector<ComputedOperator> lhs;
    BlockFunction rhs;
};

//TODO: better name
LinearSystem separate(const std::vector<ComputedOperator>& eqtn, const BCMap& bcs);

LinearSystem scale_rows(const LinearSystem& eqtn);

#endif
