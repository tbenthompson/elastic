#include "compute.h"

using namespace tbem;


MatrixOperator operator*(const MatrixOperator& op, double s) {
    std::vector<std::vector<double>> out_data(op.data.size(),
        std::vector<double>(op.data[0].size()));
    for (size_t c = 0; c < op.data.size(); c++) {
        for (size_t i = 0; i < op.data[0].size(); i++) {
            out_data[c][i] = op.data[c][i] * s;
        }
    }
    return {op.n_rows, op.n_cols, op.n_comp_rows, op.n_comp_cols, out_data};
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
std::vector<ComputedOperator>
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec& eqtn_spec)
{
    std::vector<ComputedOperator> integrals;
    for (const auto& term: eqtn_spec.terms) {
        integrals.push_back(compute_integral(bem, term));
    }
    integrals.push_back(compute_mass(bem, eqtn_spec.mass));

    return {integrals};
}

template 
std::vector<ComputedOperator>
compute_integral_equation(const BEM<2>& bem, const IntegralEquationSpec& eqtn_spec);
template 
std::vector<ComputedOperator>
compute_integral_equation(const BEM<3>& bem, const IntegralEquationSpec& eqtn_spec);

LinearSystem separate(const std::vector<ComputedOperator>& eqtn, const BCMap& bcs) {
    size_t components = eqtn[0].op.n_comp_rows;
    size_t dofs = eqtn[0].op.n_rows;
    Function rhs = constant_function(components, dofs, 0.0);
    std::vector<ComputedOperator> lhs;

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
