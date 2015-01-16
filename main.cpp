#include "spec.h"
#include "load.h"
#include "3bem/3bem.h"

using namespace tbem;

template <size_t dim>
ConstraintMatrix form_unknown_displacement_constraints(const MeshSet<dim>& meshes) 
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
QuadStrategy<dim> form_quad_strategy(const Parameters& params) {
    return QuadStrategy<dim>(
        params.obs_quad_order,
        params.src_far_quad_order,
        params.n_singular_steps,
        params.far_threshold,
        params.near_tol
    );
}

template <size_t dim>
struct BEMInput {
    const Parameters params; 
    const MeshSet<dim> meshes;
    const BCSet bcs;
    const KernelSet<dim> kernels;
    const QuadStrategy<dim> quad_strategy;
    const ConstraintMatrix unknown_displacement_constraints;

    BEMInput(const Parameters& params, const MeshSet<dim>& meshes, const BCSet& bcs):
        params(params),
        meshes(meshes),
        bcs(bcs),
        kernels(get_elastic_kernels<dim>(params.shear_modulus, params.poisson_ratio)),
        quad_strategy(form_quad_strategy<dim>(params)),
        unknown_displacement_constraints(form_unknown_displacement_constraints(meshes))
    {}
};

template <size_t dim>
BEMInput<dim> parse_into_bem_input(const std::string& filename)
{
    auto doc = parse_json(load_file(filename));
    auto params = get_parameters(doc);
    auto elements = get_elements<dim>(doc);
    auto meshes = get_meshes(elements);
    auto bcs = get_bcs(elements);

    // Setup the kernels that are necessary.
    return BEMInput<dim>(params, meshes, bcs);
}

template <size_t dim>
void compute_operator(const BEMInput<dim>& bem_input, const OperatorSpec& op_spec) {
    assert(bem_input.meshes.count(op_spec.obs_mesh) > 0);
    assert(bem_input.meshes.count(op_spec.src_mesh) > 0);
    assert(bem_input.kernels.count(op_spec.kernel) > 0);

    const auto& obs_mesh = bem_input.meshes.at(op_spec.obs_mesh);  
    const auto& src_mesh = bem_input.meshes.at(op_spec.src_mesh);
    const auto& kernel = bem_input.kernels.at(op_spec.kernel);

    auto problem = make_problem(src_mesh, obs_mesh, *kernel);
    auto op = mesh_to_mesh_operator(problem, bem_input.quad_strategy);

    std::cout << op_spec << std::endl;
    std::cout << op.rows << std::endl;
    std::cout << op.cols << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto bem_input = parse_into_bem_input<2>(filename);
    auto displacement_BIE = get_displacement_BIE();
    for (const auto& term: displacement_BIE.terms) {
        compute_operator(bem_input, term);
    }
    auto traction_BIE = get_traction_BIE();
    for (const auto& term: traction_BIE.terms) {
        compute_operator(bem_input, term);
    }
}
