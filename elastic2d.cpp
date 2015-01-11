#include "3bem/3bem.h"
#include "3bem/elastic_kernels.h"
#include "elastic.h"

using namespace tbem;

std::string remove_extension(const std::string& filename) 
{
    size_t last_dot = filename.find_last_of(".");
    if (last_dot == std::string::npos) {
        return filename;
    }
    return filename.substr(0, last_dot); 
}

ConstraintMatrix form_constraint_matrix(const Mesh<2>& surface, const Mesh<2>& fault) 
{
    auto continuity = mesh_continuity(surface.begin());
    auto cut_continuity = cut_at_intersection(
        continuity, surface.begin(), fault.begin()
    );
    auto constraints = convert_to_constraints(cut_continuity);
    auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);
    return constraint_matrix;
}

template <size_t dim>
struct KernelSet {
    const ElasticDisplacement<dim> displacement;
    const ElasticTraction<dim> traction;
    const ElasticAdjointTraction<dim> adjoint_traction;
    const ElasticHypersingular<dim> hypersingular;

    KernelSet(double shear_modulus, double poisson_ratio):
        displacement(shear_modulus, poisson_ratio),
        traction(shear_modulus, poisson_ratio),
        adjoint_traction(shear_modulus, poisson_ratio),
        hypersingular(shear_modulus, poisson_ratio)
    {}
};

template <size_t dim>
struct Solver {
    const Parameters params; 
    const MeshesAndBCs<dim> meshes_bcs;
    const KernelSet<dim> kernels;
    const QuadStrategy<dim> quad_strategy;
    const ConstraintMatrix constraint_matrix;

    Solver(const Parameters& params, const MeshesAndBCs<dim>& meshes_bcs):
        params(params),
        meshes_bcs(meshes_bcs),
        kernels(params.shear_modulus, params.poisson_ratio),
        quad_strategy(
            params.obs_quad_order,
            params.src_far_quad_order,
            params.n_singular_steps,
            params.far_threshold,
            params.near_tol
        ),
        constraint_matrix(
            form_constraint_matrix(meshes_bcs.traction_mesh, meshes_bcs.slip_mesh)
        )
    {}
};

template <size_t dim>
Solver<dim> parse_into_solver(const std::string& filename) {
    auto doc = parse_json(load_file(filename));
    auto params = get_parameters(doc);
    auto meshes_bcs = get_meshes_bcs(doc);

    // Setup the kernels that are necessary.
    return Solver<dim>(params, meshes_bcs);
}

template <size_t dim>
std::vector<Vec<double,dim>> compute_traction_rhs(const Solver<dim>& solver) {
    const auto& disp_mesh = solver.meshes_bcs.displacement_mesh;
    const auto& trac_mesh = solver.meshes_bcs.traction_mesh;
    const auto& slip_mesh = solver.meshes_bcs.slip_mesh;

    auto p_mass = make_problem(
        trac_mesh, trac_mesh, IdentityTensor<dim,dim,dim>(),
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.traction_bcs.facets)
    );
    auto rhs_mass = mass_term(p_mass, solver.quad_strategy);

    auto p_slip_trac = make_problem(
        slip_mesh, trac_mesh, solver.kernels.hypersingular,
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.slip_bcs.facets)
    );
    auto int_slip_trac = direct_interact(p_slip_trac, solver.quad_strategy);

    auto p_trac_trac = make_problem(
        trac_mesh, trac_mesh, solver.kernels.adjoint_traction,
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.traction_bcs.facets)
    );
    auto int_trac_trac = direct_interact(p_trac_trac, solver.quad_strategy);

    auto p_disp_trac = make_problem(
        disp_mesh, trac_mesh, solver.kernels.hypersingular,
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.displacement_bcs.facets)
    );
    auto int_disp_trac = direct_interact(p_disp_trac, solver.quad_strategy);
                                
    std::vector<Vec<double,dim>> trac_rhs_all_dofs(trac_mesh.n_dofs());
    for (unsigned int i = 0; i < trac_rhs_all_dofs.size(); i++) {
        trac_rhs_all_dofs[i] = int_slip_trac[i] 
            - int_trac_trac[i] 
            + int_disp_trac[i]
            - rhs_mass[i];
    }
    auto reduced_trac_rhs = solver.constraint_matrix.get_reduced(trac_rhs_all_dofs);
    return reduced_trac_rhs;
}

template <size_t dim>
std::vector<Vec<double,dim>> compute_displacement_rhs(const Solver<dim>& solver) {
    const auto& disp_mesh = solver.meshes_bcs.displacement_mesh;
    const auto& trac_mesh = solver.meshes_bcs.traction_mesh;

    auto p_mass = make_problem(
        disp_mesh, disp_mesh, IdentityTensor<dim,dim,dim>(),
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.displacement_bcs.facets)
    );
    auto rhs_mass = mass_term(p_mass, solver.quad_strategy);

    auto p_trac_disp = make_problem(
        trac_mesh, disp_mesh, solver.kernels.displacement,
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.traction_bcs.facets)
    );
    auto int_trac_disp = direct_interact(p_trac_disp, solver.quad_strategy);

    auto p_disp_disp = make_problem(
        disp_mesh, disp_mesh, solver.kernels.traction,
        reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.displacement_bcs.facets)
    );
    auto int_disp_disp = direct_interact(p_disp_disp, solver.quad_strategy);
                                
    std::vector<Vec<double,dim>> disp_rhs_all_dofs(disp_mesh.n_dofs());
    for (unsigned int i = 0; i < disp_rhs_all_dofs.size(); i++) {
        disp_rhs_all_dofs[i] = int_trac_disp[i] - int_disp_disp[i] - rhs_mass[i];
    }

    auto reduced_disp_rhs = solver.constraint_matrix.get_reduced(disp_rhs_all_dofs);
    return reduced_disp_rhs;
}

template <size_t dim>
struct RHSVector
{
    std::vector<Vec<double,dim>> trac;
    std::vector<Vec<double,dim>> disp;
};

template <size_t dim>
RHSVector<dim> compute_rhs(const Solver<dim>& solver) {
    auto traction_rhs = compute_traction_rhs(solver);
    auto displacement_rhs = compute_displacement_rhs(solver);
    return {traction_rhs, displacement_rhs};
}

template <size_t dim>
struct LHSMatrices
{
    typedef std::vector<Vec<Vec<double,dim>,dim>> MatType;
    const MatType trac_trac;
    const MatType disp_trac;
    const MatType trac_disp;
    const MatType disp_disp;
};

template <size_t dim>
LHSMatrices<dim> compute_lhs_matrices(const Solver<dim>& solver) {
    const auto& disp_mesh = solver.meshes_bcs.displacement_mesh;
    const auto& trac_mesh = solver.meshes_bcs.traction_mesh;

    auto p_trac_trac = make_problem(
        trac_mesh, trac_mesh, solver.kernels.hypersingular, {}
    );
    auto trac_trac = interact_matrix(p_trac_trac, solver.quad_strategy);
    for (size_t i = 0; i < trac_trac.size(); i++) {
        trac_trac[i] = -trac_trac[i];
    }

    auto p_disp_trac = make_problem(
        disp_mesh, trac_mesh, solver.kernels.adjoint_traction, {}
    );
    auto disp_trac = interact_matrix(p_disp_trac, solver.quad_strategy);

    auto p_trac_disp = make_problem(
        trac_mesh, disp_mesh, solver.kernels.traction, {}
    );
    auto trac_disp = interact_matrix(p_trac_disp, solver.quad_strategy);

    auto p_disp_disp = make_problem(
        disp_mesh, disp_mesh, solver.kernels.displacement, {}
    );
    auto disp_disp = interact_matrix(p_disp_disp, solver.quad_strategy);
    for (size_t i = 0; i < disp_disp.size(); i++) {
        disp_disp[i] = -disp_disp[i];
    }

    return {trac_trac, disp_trac, trac_disp, disp_disp};
}

std::vector<double> concatenate(const std::vector<double>& A,
                                const std::vector<double>& B)
{
    std::vector<double> AB;
    AB.reserve(A.size() + B.size());
    AB.insert(AB.end(), A.begin(), A.end());
    AB.insert(AB.end(), B.begin(), B.end()); 
    return AB;
}

template <size_t dim>
std::vector<Vec<double,dim>> solve_bem_linear_system(const Solver<dim>& solver,
    const LHSMatrices<dim>& lhs, const RHSVector<dim>& rhs)
{
    int trac_dofs = solver.meshes_bcs.traction_mesh.n_dofs();
    int disp_dofs = solver.meshes_bcs.displacement_mesh.n_dofs();
            
    int reduced_trac_dofs = rhs.trac.size();

    auto rhs_all_meshes = concatenate(
        reinterpret_vector<double>(rhs.trac),
        reinterpret_vector<double>(rhs.disp)
    );

    int count = 0;
    auto reduced_soln = solve_system(rhs_all_meshes, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            count++;
            auto vec_x = reinterpret_vector<Vec<double,dim>>(x);
            auto trac_begin = vec_x.begin();
            auto disp_begin = vec_x.begin() + reduced_trac_dofs;
            std::vector<Vec<double,dim>> reduced_trac(trac_begin, disp_begin);
            std::vector<Vec<double,dim>> reduced_disp(disp_begin, vec_x.end());
            auto all_trac = solver.constraint_matrix.get_all(reduced_trac, trac_dofs);
            auto all_disp = solver.constraint_matrix.get_all(reduced_disp, disp_dofs);

            auto eval_trac_trac = bem_mat_mult(
                lhs.trac_trac, solver.kernels.hypersingular, trac_dofs, all_trac
            );
            auto eval_disp_trac = bem_mat_mult(
                lhs.disp_trac, solver.kernels.adjoint_traction, trac_dofs, all_disp
            );
            std::vector<Vec<double,dim>> eval_trac(trac_dofs);
            for (size_t i = 0; i < trac_dofs; i++) {
                eval_trac[i] = eval_trac_trac[i] + eval_disp_trac[i];
            }
            auto reduced_eval_trac = solver.constraint_matrix.get_reduced(eval_trac);

            auto eval_trac_disp = bem_mat_mult(
                lhs.trac_disp, solver.kernels.traction, disp_dofs, all_trac
            );
            auto eval_disp_disp = bem_mat_mult(
                lhs.disp_disp, solver.kernels.displacement, disp_dofs, all_disp
            );
            std::vector<Vec<double,dim>> eval_disp(disp_dofs);
            for (size_t i = 0; i < disp_dofs; i++) {
                eval_disp[i] = eval_trac_disp[i] + eval_disp_disp[i];
            }
            auto reduced_eval_disp = solver.constraint_matrix.get_reduced(eval_disp);

            auto result = concatenate(
                reinterpret_vector<double>(reduced_eval_trac),
                reinterpret_vector<double>(reduced_eval_disp)
            );
            std::copy(result.begin(), result.end(), y.begin());
        }
    );
    std::cout << "Took " << count << " iterations." << std::endl;
    return reinterpret_vector<Vec<double,dim>>(reduced_soln);
}

template <size_t dim>
void output_solution(const Solver<dim>& solver,
    const std::vector<Vec<double,dim>>& reduced_soln, 
    int reduced_trac_dofs,
    const std::string& in_filename)
{
    int trac_dofs = solver.meshes_bcs.traction_mesh.n_dofs();
    int disp_dofs = solver.meshes_bcs.displacement_mesh.n_dofs();

    auto trac_begin = reduced_soln.begin();
    auto disp_begin = reduced_soln.begin() + reduced_trac_dofs;
    std::vector<Vec<double,dim>> reduced_trac(trac_begin, disp_begin);
    std::vector<Vec<double,dim>> reduced_disp(disp_begin, reduced_soln.end());
    auto all_trac_soln = solver.constraint_matrix.get_all(reduced_trac, trac_dofs);
    auto all_disp_soln = solver.constraint_matrix.get_all(reduced_disp, disp_dofs);

    auto in_filename_root = remove_extension(in_filename);
    auto out_filename_trac = in_filename_root + ".trac_out";
    auto out_filename_disp = in_filename_root + ".disp_out";
    if (trac_dofs > 0) {
        auto file = HDFOutputter(in_filename_root + ".trac_out");
        out_surface(file, solver.meshes_bcs.traction_mesh, all_trac_soln, dim);
        std::cout << all_trac_soln.size() << std::endl;
    }

    if (disp_dofs > 0) {
        auto file = HDFOutputter(in_filename_root + ".disp_out");
        out_surface(file, solver.meshes_bcs.displacement_mesh, all_disp_soln, dim);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto solver = parse_into_solver<2>(filename);

    std::cout << "Number of displacement DOFs: " 
        << solver.meshes_bcs.displacement_mesh.n_dofs() << std::endl;
    std::cout << "Number of traction DOFs: " 
        << solver.meshes_bcs.traction_mesh.n_dofs() << std::endl;
    std::cout << "Number of slip DOFs: " 
        << solver.meshes_bcs.slip_mesh.n_dofs() << std::endl;

    auto rhs = compute_rhs(solver);
    auto lhs = compute_lhs_matrices(solver);
    auto soln = solve_bem_linear_system(solver, lhs, rhs);
    output_solution(solver, soln, rhs.trac.size(), filename);
}
