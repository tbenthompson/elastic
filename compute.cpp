#include "compute.h"

using namespace tbem;

template <size_t dim>
BEM<dim> parse_into_bem(const std::string& filename)
{
    auto doc = parse_json(load_file(filename));
    auto params = get_parameters(doc);
    auto elements = get_elements<dim>(doc);
    auto meshes = get_meshes(elements);
    auto bcs = get_bcs(elements);

    // Setup the kernels that are necessary.
    return BEM<dim>(params, meshes, bcs);
}

template 
BEM<2> parse_into_bem(const std::string& filename);
template 
BEM<3> parse_into_bem(const std::string& filename);

template <size_t dim>
BEM<dim>::BEM(const Parameters& params,
        const MeshMap<dim>& meshes, const BCMap& bcs):
    meshes(meshes),
    bcs(bcs),
    kernels(get_elastic_kernels<dim>(params.shear_modulus, params.poisson_ratio)),
    quad_strategy(form_quad_strategy(params)),
    displacement_constraints(form_displacement_constraints(meshes))
{}

template <size_t dim>
ConstraintMatrix BEM<dim>::form_displacement_constraints(const MeshMap<dim>& meshes) 
{
    auto continuity = mesh_continuity(meshes.at("traction").begin());
    auto cut_continuity = cut_at_intersection(
        continuity, meshes.at("traction").begin(), meshes.at("slip").begin()
    );
    auto constraints = convert_to_constraints(cut_continuity);
    auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);
    return constraint_matrix;
}

template <size_t dim>
QuadStrategy<dim> BEM<dim>::form_quad_strategy(const Parameters& params) {
    return QuadStrategy<dim>(
        params.obs_quad_order,
        params.src_far_quad_order,
        params.n_singular_steps,
        params.far_threshold,
        params.near_tol
    );
}

MatrixOperator operator*(const MatrixOperator& op, double s) {
    std::vector<std::vector<double>> out_data(op.data.size(),
        std::vector<double>(op.data[0].size()));
    for (size_t c = 0; c < op.data.size(); c++) {
        for (size_t i = 0; i < op.data[0].size(); i++) {
            out_data[c][i] = op.data[c][i] * s;
        }
    }
    return {op.rows, op.cols, op.n_comp_rows, op.n_comp_cols, out_data};
}

template <size_t dim>
ComputedOperator 
BEM<dim>::compute_integral(const IntegralSpec& op_spec) {
    assert(meshes.count(op_spec.obs_mesh) > 0);
    assert(meshes.count(op_spec.src_mesh) > 0);
    assert(kernels.count(op_spec.kernel) > 0);

    const auto& obs_mesh = meshes.at(op_spec.obs_mesh);  
    const auto& src_mesh = meshes.at(op_spec.src_mesh);
    const auto& kernel = kernels.at(op_spec.kernel);

    auto problem = make_problem(src_mesh, obs_mesh, *kernel);
    auto op = mesh_to_mesh_operator(problem, quad_strategy);

    return {op * op_spec.multiplier, op_spec.src_mesh, op_spec.function};
}

template <size_t dim>
ComputedOperator
BEM<dim>::compute_mass(const MassSpec& op_spec) {
    assert(meshes.count(op_spec.obs_mesh) > 0);

    const auto& obs_mesh = meshes.at(op_spec.obs_mesh);  

    IdentityTensor<dim,dim,dim> identity;
    auto problem = make_problem(obs_mesh, obs_mesh, identity);
    auto mass_op = mass_operator(problem, quad_strategy);

    return {mass_op * op_spec.multiplier, op_spec.obs_mesh, op_spec.function};
}


template <size_t dim>
ComputedIntegralEquation
BEM<dim>::compute_integral_equation(const IntegralEquationSpec& eqtn_spec) 
{
    std::vector<ComputedOperator> integrals;
    for (const auto& term: eqtn_spec.terms) {
        integrals.push_back(compute_integral(term));
    }
    integrals.push_back(compute_mass(eqtn_spec.mass));

    return {integrals};
}

LinearSystem separate(const ComputedIntegralEquation& eqtn, const BCMap& bcs) {
    size_t components = eqtn.terms[0].op.n_comp_rows;
    size_t dofs = eqtn.terms[0].op.rows;
    Function rhs = constant_function(components, dofs, 0.0);
    std::vector<ComputedOperator> lhs;

    for (const auto& term: eqtn.terms) {
        FieldDescriptor field_desc{term.src_mesh, term.function};
        auto it = bcs.find(field_desc);
        if (it == bcs.end()) {
            lhs.push_back(term);
        } else 
        {
            const auto& bc = it->second;
            //negate because the term is shifted to the other side of the equation.
            std::cout << term.src_mesh << term.function << std::endl;
            rhs -= apply_operator(term.op, bc);
        }
    }

    return {lhs, rhs};
}

MatrixOperator divide_rows(const MatrixOperator& op, const Function& rhs) {
    std::vector<std::vector<double>> out_data(op.data.size(),
            std::vector<double>(op.data[0].size()));
    for (size_t d0 = 0; d0 < op.n_comp_rows; d0++) {
        for (size_t d1 = 0; d1 < op.n_comp_cols; d1++) {
            for (size_t i = 0; i < op.rows; i++) {
                for (size_t j = 0; j < op.cols; j++) {
                    auto val = 
                        op.data[d0 * op.n_comp_cols + d1][i * op.cols + j] / rhs[d0][i];
                    out_data[d0 * op.n_comp_cols + d1][i * op.cols + j] = val;
                }
            }
        }
    }
    return {op.rows, op.cols, op.n_comp_rows, op.n_comp_cols, out_data};
}

LinearSystem scale_rows(const LinearSystem& eqtn) 
{
    std::vector<ComputedOperator> integrals;
    for (const auto& term: eqtn.lhs) {
        integrals.push_back({
            divide_rows(term.op, eqtn.rhs),
            term.src_mesh,
            term.function
        });
    }
    return {
        integrals,
        constant_function(eqtn.rhs.size(), eqtn.rhs[0].size(), 1.0)
    };
}

template struct BEM<2>;
template struct BEM<3>;
