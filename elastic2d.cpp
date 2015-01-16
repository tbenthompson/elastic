#include "3bem/3bem.h"
#include "3bem/elastic_kernels.h"
#include "load.h"

using namespace tbem;

// std::string remove_extension(const std::string& filename) 
// {
//     size_t last_dot = filename.find_last_of(".");
//     if (last_dot == std::string::npos) {
//         return filename;
//     }
//     return filename.substr(0, last_dot); 
// }
// 
// {
//     auto continuity = mesh_continuity(surface.begin());
//     auto cut_continuity = cut_at_intersection(
//         continuity, surface.begin(), fault.begin()
//     );
//     auto constraints = convert_to_constraints(cut_continuity);
//     auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);
//     return constraint_matrix;
// }
// 
// template <size_t dim>
// struct KernelSet {
//     const ElasticDisplacement<dim> displacement;
//     const ElasticTraction<dim> traction;
//     const ElasticAdjointTraction<dim> adjoint_traction;
//     const ElasticHypersingular<dim> hypersingular;
// 
//     KernelSet(double shear_modulus, double poisson_ratio):
//         displacement(shear_modulus, poisson_ratio),
//         traction(shear_modulus, poisson_ratio),
//         adjoint_traction(shear_modulus, poisson_ratio),
//         hypersingular(shear_modulus, poisson_ratio)
//     {}
// };
// 
// template <size_t dim>
// struct Solver {
//     const Parameters params; 
//     const MeshesAndBCs<dim> meshes_bcs;
//     const KernelSet<dim> kernels;
//     const QuadStrategy<dim> quad_strategy;
//     const ConstraintMatrix constraint_matrix;
// 
//     Solver(const Parameters& params, const MeshesAndBCs<dim>& meshes_bcs):
//         params(params),
//         meshes_bcs(meshes_bcs),
//         kernels(params.shear_modulus, params.poisson_ratio),
//         quad_strategy(
//             params.obs_quad_order,
//             params.src_far_quad_order,
//             params.n_singular_steps,
//             params.far_threshold,
//             params.near_tol
//         ),
//         constraint_matrix(
//             form_constraint_matrix(meshes_bcs.traction_mesh, meshes_bcs.slip_mesh)
//         )
//     {}
// };
// 
// template <size_t dim>
// Solver<dim> parse_into_solver(const std::string& filename) {
//     auto doc = parse_json(load_file(filename));
//     auto params = get_parameters(doc);
//     auto meshes_bcs = get_meshes_bcs(doc);
// 
//     // Setup the kernels that are necessary.
//     return Solver<dim>(params, meshes_bcs);
// }
// 
// 
// template <size_t dim>
// std::vector<Vec<double,dim>> compute_traction_rhs(const Solver<dim>& solver) {
//     const auto& disp_mesh = solver.meshes_bcs.displacement_mesh;
//     const auto& trac_mesh = solver.meshes_bcs.traction_mesh;
//     const auto& slip_mesh = solver.meshes_bcs.slip_mesh;
// 
//     auto p_mass = make_problem(
//         trac_mesh, trac_mesh, IdentityTensor<dim,dim,dim>(),
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.traction_bcs.facets)
//     );
//     auto rhs_mass = mass_term(p_mass, solver.quad_strategy);
// 
//     auto p_slip_trac = make_problem(
//         slip_mesh, trac_mesh, solver.kernels.hypersingular,
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.slip_bcs.facets)
//     );
//     auto int_slip_trac = direct_interact(p_slip_trac, solver.quad_strategy);
// 
//     auto p_trac_trac = make_problem(
//         trac_mesh, trac_mesh, solver.kernels.adjoint_traction,
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.traction_bcs.facets)
//     );
//     auto int_trac_trac = direct_interact(p_trac_trac, solver.quad_strategy);
// 
//     auto p_disp_trac = make_problem(
//         disp_mesh, trac_mesh, solver.kernels.hypersingular,
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.displacement_bcs.facets)
//     );
//     auto int_disp_trac = direct_interact(p_disp_trac, solver.quad_strategy);
//                                 
//     std::vector<Vec<double,dim>> trac_rhs_all_dofs(trac_mesh.n_dofs());
//     for (unsigned int i = 0; i < trac_rhs_all_dofs.size(); i++) {
//         trac_rhs_all_dofs[i] = int_slip_trac[i] 
//             - int_trac_trac[i] 
//             + int_disp_trac[i]
//             - rhs_mass[i];
//     }
//     auto reduced_trac_rhs = solver.constraint_matrix.get_reduced(trac_rhs_all_dofs);
//     return reduced_trac_rhs;
// }
// 
// template <size_t dim>
// std::vector<Vec<double,dim>> compute_displacement_rhs(const Solver<dim>& solver) {
//     const auto& disp_mesh = solver.meshes_bcs.displacement_mesh;
//     const auto& trac_mesh = solver.meshes_bcs.traction_mesh;
// 
//     auto p_mass = make_problem(
//         disp_mesh, disp_mesh, IdentityTensor<dim,dim,dim>(),
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.displacement_bcs.facets)
//     );
//     auto rhs_mass = mass_term(p_mass, solver.quad_strategy);
// 
//     auto p_trac_disp = make_problem(
//         trac_mesh, disp_mesh, solver.kernels.displacement,
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.traction_bcs.facets)
//     );
//     auto int_trac_disp = direct_interact(p_trac_disp, solver.quad_strategy);
// 
//     auto p_disp_disp = make_problem(
//         disp_mesh, disp_mesh, solver.kernels.traction,
//         reinterpret_vector<Vec<double,dim>>(solver.meshes_bcs.displacement_bcs.facets)
//     );
//     auto int_disp_disp = direct_interact(p_disp_disp, solver.quad_strategy);
//                                 
//     std::vector<Vec<double,dim>> disp_rhs_all_dofs(disp_mesh.n_dofs());
//     for (unsigned int i = 0; i < disp_rhs_all_dofs.size(); i++) {
//         disp_rhs_all_dofs[i] = int_trac_disp[i] - int_disp_disp[i] - rhs_mass[i];
//     }
// 
//     return disp_rhs_all_dofs;
// }
// 
// template <size_t dim>
// struct RHSVector
// {
//     std::vector<Vec<double,dim>> trac;
//     std::vector<Vec<double,dim>> disp;
// };
// 
// template <size_t dim>
// RHSVector<dim> compute_rhs(const Solver<dim>& solver) {
//     auto traction_rhs = compute_traction_rhs(solver);
//     auto displacement_rhs = compute_displacement_rhs(solver);
//     return {traction_rhs, displacement_rhs};
// }
// 
// template <size_t dim>
// struct LHSMatrices
// {
//     typedef std::vector<Vec<Vec<double,dim>,dim>> MatType;
//     const MatType disp_disp;
//     const MatType trac_disp;
//     const MatType disp_trac;
//     const MatType trac_trac;
// };
// 
// template <size_t dim>
// LHSMatrices<dim> compute_lhs_matrices(const Solver<dim>& solver) {
//     const auto& unknown_trac_mesh = solver.meshes_bcs.displacement_mesh;
//     const auto& unknown_disp_mesh = solver.meshes_bcs.traction_mesh;
// 
//     auto p_disp_disp = make_problem(
//         unknown_disp_mesh, unknown_disp_mesh, solver.kernels.hypersingular, {}
//     );
//     auto disp_disp = interact_matrix(p_disp_disp, solver.quad_strategy);
//     for (size_t i = 0; i < disp_disp.size(); i++) {
//         disp_disp[i] = -disp_disp[i];
//     }
// 
//     auto p_trac_disp = make_problem(
//         unknown_trac_mesh, unknown_disp_mesh, solver.kernels.adjoint_traction, {}
//     );
//     auto trac_disp = interact_matrix(p_trac_disp, solver.quad_strategy);
// 
//     auto p_disp_trac = make_problem(
//         unknown_disp_mesh, unknown_trac_mesh, solver.kernels.traction, {}
//     );
//     auto disp_trac = interact_matrix(p_trac_disp, solver.quad_strategy);
// 
//     auto p_trac_trac = make_problem(
//         unknown_trac_mesh, unknown_trac_mesh, solver.kernels.displacement, {}
//     );
//     auto trac_trac = interact_matrix(p_trac_trac, solver.quad_strategy);
//     for (size_t i = 0; i < trac_trac.size(); i++) {
//         trac_trac[i] = -trac_trac[i];
//     }
// 
//     return {disp_disp, trac_disp, disp_trac, trac_trac};
// }
// 
// std::vector<double> concatenate(const std::vector<double>& A,
//                                 const std::vector<double>& B)
// {
//     std::vector<double> AB;
//     AB.reserve(A.size() + B.size());
//     AB.insert(AB.end(), A.begin(), A.end());
//     AB.insert(AB.end(), B.begin(), B.end()); 
//     return AB;
// }
// 
// template <size_t dim>
// RHSVector<dim> extract_and_expand_fields(
//     const Solver<dim>& solver,
//     const std::vector<Vec<double,dim>> all_dofs,
//     int reduced_trac_dofs)
// {
//     int trac_dofs = solver.meshes_bcs.traction_mesh.n_dofs();
// 
//     auto trac_begin = all_dofs.begin();
//     auto disp_begin = all_dofs.begin() + reduced_trac_dofs;
//     std::vector<Vec<double,dim>> reduced_trac(trac_begin, disp_begin);
//     std::vector<Vec<double,dim>> all_disp(disp_begin, all_dofs.end());
// 
//     return {
//         solver.constraint_matrix.get_all(reduced_trac, trac_dofs),
//         all_disp
//     };
// }
// 
// template <size_t dim>
// RHSVector<dim> solve_bem_linear_system(const Solver<dim>& solver,
//     const LHSMatrices<dim>& lhs, const RHSVector<dim>& rhs)
// {
//     int trac_dofs = solver.meshes_bcs.traction_mesh.n_dofs();
//     int disp_dofs = solver.meshes_bcs.displacement_mesh.n_dofs();
//             
//     int reduced_trac_dofs = rhs.trac.size();
// 
//     auto rhs_all_meshes = concatenate(
//         reinterpret_vector<double>(rhs.trac),
//         reinterpret_vector<double>(rhs.disp)
//     );
// 
//     int count = 0;
//     auto reduced_soln = solve_system(rhs_all_meshes, 1e-5,
//         [&] (std::vector<double>& x, std::vector<double>& y) {
//             count++;
//             auto vec_x = reinterpret_vector<Vec<double,dim>>(x);
//             auto fields = extract_and_expand_fields(solver, vec_x, reduced_trac_dofs);
//             const auto& unknown_trac = fields.disp;
//             const auto& unknown_disp = fields.trac;
// 
//             auto eval_trac_trac = bem_mat_mult(
//                 lhs.trac_trac, solver.kernels.hypersingular,
//                 disp_dofs, unknown_trac
//             );
//             auto eval_disp_trac = bem_mat_mult(
//                 lhs.disp_trac, solver.kernels.adjoint_traction,
//                 disp_dofs, unknown_disp
//             );
//             std::vector<Vec<double,dim>> eval_trac(disp_dofs);
//             for (size_t i = 0; i < disp_dofs; i++) {
//                 eval_trac[i] = eval_trac_trac[i] + eval_disp_trac[i];
//             }
// 
//             auto eval_trac_disp = bem_mat_mult(
//                 lhs.trac_disp, solver.kernels.traction,
//                 trac_dofs, unknown_trac
//             );
//             auto eval_disp_disp = bem_mat_mult(
//                 lhs.disp_disp, solver.kernels.displacement,
//                 trac_dofs, unknown_disp
//             );
//             std::vector<Vec<double,dim>> eval_disp(trac_dofs);
//             for (size_t i = 0; i < trac_dofs; i++) {
//                 eval_disp[i] = eval_trac_disp[i] + eval_disp_disp[i];
//             }
//             auto reduced_eval_disp = solver.constraint_matrix.get_reduced(eval_disp);
// 
//             auto result = concatenate(
//                 reinterpret_vector<double>(reduced_eval_disp),
//                 reinterpret_vector<double>(eval_trac)
//             );
//             assert(result.size() == rhs_all_meshes.size());
//             std::copy(result.begin(), result.end(), y.begin());
//         }
//     );
//     std::cout << "Took " << count << " iterations." << std::endl;
//     return extract_and_expand_fields(
//         solver, 
//         reinterpret_vector<Vec<double,dim>>(reduced_soln),
//         rhs.trac.size()
//     );
// }
// 
// template <size_t dim>
// Mesh<dim> merge_meshes(const Solver<dim>& solver) {
//     return Mesh<2>::form_union({
//         solver.meshes_bcs.displacement_mesh,
//         solver.meshes_bcs.traction_mesh
//     });
// }
// 
// template <size_t dim>
// RHSVector<dim> merge_fields(const Solver<dim>& solver, const RHSVector<dim>& soln) {
//     auto combined_disp_field = Mesh<2>::form_union({
//         solver.meshes_bcs.displacement_bcs,
//         Mesh<2>{reinterpret_vector<Facet<2>>(soln.trac)}
//     });
//     auto combined_trac_field = Mesh<2>::form_union({
//         Mesh<2>{reinterpret_vector<Facet<2>>(soln.disp)},
//         solver.meshes_bcs.traction_bcs
//     });
//     return {
//         reinterpret_vector<Vec<double,dim>>(combined_trac_field.facets),
//         reinterpret_vector<Vec<double,dim>>(combined_disp_field.facets),
//     };
// }
// 
// template <size_t dim>
// void output_solution(const Mesh<dim>& merged_mesh,
//     const RHSVector<dim>& soln, 
//     const std::string& in_filename)
// {
//     auto in_filename_root = remove_extension(in_filename);
//     auto out_filename_trac = in_filename_root + ".trac_out";
//     auto out_filename_disp = in_filename_root + ".disp_out";
//     out_surface(HDFOutputter(out_filename_trac), merged_mesh, soln.trac, dim);
//     out_surface(HDFOutputter(out_filename_disp), merged_mesh, soln.disp, dim);
// }
// 
// 
// template <size_t dim>
// std::vector<Vec<double,dim>> interior_eval(const Solver<dim>& solver,
//     const Mesh<dim>& mesh, const RHSVector<dim>& fields,
//     const std::vector<Vec<double,dim>>& locs) 
// {
//     //u(x) = (field=traction, K=displacement) - 
//     //       (field=displacement, K=traction)
//     //t(x) = (field=traction, K=adjoint_traction) - 
//     //       (field=displacement, K=hypersingular)
// 
//     std::vector<Vec<double,dim>> results(locs.size());
// #pragma omp parallel for
//     for (size_t i = 0; i < locs.size(); i++) {
//         auto pt = locs[i];
//         auto dir_to_origin = -pt;
//         ObsPt<dim> obs_pt = {0.001, pt, dir_to_origin, dir_to_origin};
// 
//         auto p_trac = make_problem(
//             mesh,
//             Mesh<dim>{{}},
//             solver.kernels.displacement,
//             fields.trac
//         );
//         auto eval_trac = eval_integral_equation(
//             p_trac,
//             solver.quad_strategy,
//             obs_pt
//         );
// 
//         auto p_disp = make_problem(
//             mesh,
//             Mesh<dim>{{}},
//             solver.kernels.traction,
//             fields.disp
//         );
//         auto eval_disp = eval_integral_equation(
//             p_disp,
//             solver.quad_strategy,
//             obs_pt
//         );
//         results[i] = eval_disp - eval_trac;
//     }
//     return results;
// }
// 
// template <size_t dim>
// void output_interior(const std::vector<Vec<double,dim>>& pts,
//     const std::vector<Vec<double,dim>>& values,
//     const std::string& in_filename)
// {
//     auto in_filename_root = remove_extension(in_filename);
//     auto out_filename_disp = in_filename_root + ".interior_disp_out";
//     out_volume(HDFOutputter(out_filename_disp), pts, values, dim);
// }
// 
int main(int argc, char* argv[]) {
    // if (argc < 2) {
    //     std::cout << "Usage is 'elastic_process filename'" << std::endl;
    //     return 1;
    // }
    // 
    // auto filename = argv[1];
    // auto solver = parse_into_solver<2>(filename);

    // std::cout << "Number of displacement DOFs: " 
    //     << solver.meshes_bcs.displacement_mesh.n_dofs() << std::endl;
    // std::cout << "Number of traction DOFs: " 
    //     << solver.meshes_bcs.traction_mesh.n_dofs() << std::endl;
    // std::cout << "Number of slip DOFs: " 
    //     << solver.meshes_bcs.slip_mesh.n_dofs() << std::endl;

    // TIC;
    // auto rhs = compute_rhs(solver);
    // auto lhs = compute_lhs_matrices(solver);
    // TOC("Forming the linear system");
    // TIC2;
    // auto soln = solve_bem_linear_system(solver, lhs, rhs);
    // TOC("Solve the linear system");
    // TIC2;
    // auto merged_mesh = merge_meshes(solver);
    // auto merged_fields = merge_fields(solver, soln);
    // output_solution(merged_mesh, merged_fields, filename);
    // TOC("Outputting the solution");

    // TIC2;
    // auto x_vals = linspace(-1, 1, 20);
    // auto y_vals = linspace(-1, 1, 20);
    // std::vector<Vec<double,2>> locs;
    // for (size_t i = 0; i < x_vals.size(); i++) {
    //     std::cout << x_vals[i] << std::endl;
    //     for (size_t j = 0; j < y_vals.size(); j++) {
    //         locs.push_back({{x_vals[i], y_vals[j]}});
    //     }
    // }

    // auto interior_disp = interior_eval(solver, merged_mesh, merged_fields, locs);
    // output_interior(locs, interior_disp, filename);
    // TOC("Evaluating at " + std::to_string(locs.size()) + " interior points");
}
