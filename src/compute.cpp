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

    return {
        op * op_spec.multiplier,
        op_spec.obs_mesh,
        op_spec.src_mesh,
        op_spec.function
    };
}

template <size_t dim>
ComputedOperator compute_mass(const BEM<dim>& bem, const MassSpec& op_spec) {
    assert(bem.meshes.count(op_spec.obs_mesh) > 0);

    const auto& obs_mesh = bem.meshes.at(op_spec.obs_mesh);  

    IdentityTensor<dim,dim,dim> identity;
    auto problem = make_problem(obs_mesh, obs_mesh, identity);
    auto mass_op = mass_operator(problem, bem.quad_strategy);

    return {
        mass_op * op_spec.multiplier,
        op_spec.obs_mesh,
        op_spec.obs_mesh,
        op_spec.function
    };
}

template <size_t dim>
void print_operator(std::ostream& os, const BEM<dim>& bem_input, const IntegralSpec& term) 
{
    if (bem_input.meshes.at(term.obs_mesh).n_facets() == 0) {
        return;
    }
    if (bem_input.meshes.at(term.src_mesh).n_facets() == 0) {
        return;
    }
    os << "Computing operator: " << std::endl;
    os << "    Observation mesh: " << term.obs_mesh << std::endl;
    os << "    Source mesh: " << term.src_mesh << std::endl;
    os << "    Kernel: " << term.kernel << std::endl;
    os << "    Function: " << term.function << std::endl; 
}

template <size_t dim>
ComputedEquation
compute_integral_equation(const BEM<dim>& bem, const IntegralEquationSpec<dim>& eqtn_spec)
{
    ComputedEquation integrals;
    for (const auto& term: eqtn_spec.terms) {
        print_operator(std::cout, bem, term);
        integrals.push_back(compute_integral(bem, term));
    }
    integrals.push_back(compute_mass(bem, eqtn_spec.mass));

    return {integrals};
}

template 
ComputedEquation
compute_integral_equation(const BEM<2>& bem, const IntegralEquationSpec<2>& eqtn_spec);
template 
ComputedEquation
compute_integral_equation(const BEM<3>& bem, const IntegralEquationSpec<3>& eqtn_spec);

LinearSystem 
evaluate_computable_terms(const ComputedEquation& eqtn, const FunctionMap& fields)
{
    size_t components = eqtn[0].op.n_comp_rows;
    size_t dofs = eqtn[0].op.ops[0].n_rows;
    BlockFunction evaluated_terms = constant_function(components, dofs, 0.0);
    ComputedEquation lhs;

    for (const auto& term: eqtn) {
        FieldDescriptor field_desc{term.src_mesh, term.function};
        auto it = fields.find(field_desc);
        if (it == fields.end()) {
            lhs.push_back(term);
        } else 
        {
            const auto& bc = it->second;

            //negate because the term is shifted to the other side of the equation.
            evaluated_terms -= apply_operator(term.op, bc);
        }
    }
    return {lhs, evaluated_terms};
}
