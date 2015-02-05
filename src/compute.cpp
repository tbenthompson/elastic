#include "compute.h"
#include "spec.h"
#include "data.h"
#include "3bem/3bem.h"

using namespace tbem;

Operator operator*(const Operator& op, double s) {
    auto out = op;
    for (size_t i = 0; i < op.data.size(); i++) {
        out.data[i] *= s;
    }
    return out;
}

BlockOperator operator*(const BlockOperator& op, double s) {
    std::vector<Operator> out_data;
    for (size_t c = 0; c < op.ops.size(); c++) {
        out_data.push_back(op.ops[c] * s);
    }
    return {op.n_comp_rows, op.n_comp_cols, out_data};
}

template <size_t dim>
ComputedOperator compute_integral(const BEM<dim>& bem, const IntegralSpec& op_spec)
{
    assert(bem.meshes.count(op_spec.obs_mesh) > 0);
    assert(bem.meshes.count(op_spec.src_mesh) > 0);
    assert(bem.kernels.count(op_spec.kernel) > 0);

    const auto& obs_mesh = bem.meshes.at(op_spec.obs_mesh);  
    const auto& src_mesh = bem.meshes.at(op_spec.src_mesh);
    const auto& kernel = bem.kernels.at(op_spec.kernel);

    auto problem = make_problem(src_mesh, obs_mesh, *kernel);
    auto op = mesh_to_mesh_operator(problem, bem.quad_strategy);

    return {op * op_spec.multiplier, op_spec.src_mesh, op_spec.function};
}

template <size_t dim>
ComputedOperator compute_mass(const BEM<dim>& bem, const MassSpec& op_spec) {
    assert(bem.meshes.count(op_spec.obs_mesh) > 0);

    const auto& obs_mesh = bem.meshes.at(op_spec.obs_mesh);  

    IdentityTensor<dim,dim,dim> identity;
    auto problem = make_problem(obs_mesh, obs_mesh, identity);
    auto mass_op = mass_operator(problem, bem.quad_strategy);

    return {mass_op * op_spec.multiplier, op_spec.obs_mesh, op_spec.function};
}


template <size_t dim>
ComputedIntegralEquation
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec<dim>& eqtn_spec)
{
    ComputedIntegralEquation integrals;
    for (const auto& term: eqtn_spec.terms) {
        integrals.push_back(compute_integral(bem, term));
    }
    integrals.push_back(compute_mass(bem, eqtn_spec.mass));

    return {integrals};
}

template 
ComputedIntegralEquation
compute_integral_equation(const BEM<2>& bem, const IntegralEquationSpec<2>& eqtn_spec);
template 
ComputedIntegralEquation
compute_integral_equation(const BEM<3>& bem, const IntegralEquationSpec<3>& eqtn_spec);

LinearSystem separate(const ComputedIntegralEquation& eqtn, const BCMap& bcs) {
    size_t components = eqtn[0].op.n_comp_rows;
    size_t dofs = eqtn[0].op.ops[0].n_rows;
    BlockFunction rhs = constant_function(components, dofs, 0.0);
    ComputedIntegralEquation lhs;

    for (const auto& term: eqtn) {
        FieldDescriptor field_desc{term.src_mesh, term.function};
        auto it = bcs.find(field_desc);
        if (it == bcs.end()) {
            lhs.push_back(term);
        } else 
        {
            const auto& bc = it->second;

            //negate because the term is shifted to the other side of the equation.
            rhs -= apply_operator(term.op, bc);
        }
    }
    return {lhs, rhs};
}
